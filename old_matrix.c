#ifdef __APPLE__
#include <arm_neon.h>
#include <assert.h>
#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

typedef struct timespec Time;

Time TimeGet()
{
	Time time;
	clock_gettime(CLOCK_REALTIME, &time);
	return time;
}

double TimeDiffInSeconds(Time start, Time end)
{
	long seconds = end.tv_sec - start.tv_sec;
	long nanoseconds = end.tv_nsec - start.tv_nsec;
	return seconds + nanoseconds * 1e-9;
}

typedef struct Matrixf64
{
	size_t NRows;
	size_t NCols;
	double *Data;
} Matrixf64;

/*
Originally templated off of Savine's intro to Modern Computational Finance.
*/
Matrixf64 Matrixf64MultiplyFast(const Matrixf64 *A, const Matrixf64 *B)
{
	assert(A->NCols == B->NRows);
	Matrixf64 Result;
	Result.NRows = A->NRows;
	Result.NCols = B->NCols;
	Result.Data = calloc(Result.NRows * Result.NCols, sizeof(double));
	for (size_t i = 0; i < A->NRows; i++)
	{
		const double *ARow = A->Data + i * A->NCols;
		double *ResultRow = Result.Data + i * Result.NCols;
		for (size_t k = 0; k < B->NRows; k++)
		{
			const double *BRow = B->Data + k * B->NCols;
			const double Aik = ARow[k];
			// #pragma clang loop vectorize_width(8) // interleave_count(4)
			for (size_t j = 0; j < B->NCols; j++)
			{
				ResultRow[j] += Aik * BRow[j];
			}
		}
	}
	return Result;
}

/*
Let's write this in NEON

v1 - First verdict: compiled in debug mode, this ends up about halfway between naive mode and fast mode.
Memory access still seems to be the main bottleneck. For the next version we should actually
make this tiled.

*/
Matrixf64 Matrixf64MultiplyFaster(const Matrixf64 *A, const Matrixf64 *B)
{
	assert(A->NCols == B->NRows);
	Matrixf64 Result;
	Result.NRows = A->NRows;
	Result.NCols = B->NCols;
	Result.Data = calloc(Result.NRows * Result.NCols, sizeof(double));
	const size_t BlockSize = 8;
	const size_t NBlocks = A->NCols / BlockSize;
	const size_t Remainder = A->NCols % BlockSize;
	for (size_t i = 0; i < A->NRows; i++)
	{
		const double *ARow = A->Data + i * A->NCols;
		double *ResultRow = Result.Data + i * Result.NCols;
		for (size_t k = 0; k < B->NRows; k++)
		{
			const double *BRow = B->Data + k * B->NCols;
			// let's assume we can do 4x float64 at a time (we should be able to, NEON says 128x4 FMA)
			for (size_t j = 0; j < NBlocks; j++)
			{
				// load 4 blocks of 2 elements from Result into memory
				float64x2x4_t ResBlocks = vld4q_f64(ResultRow + j * BlockSize);
				// load 4 blocks of 2 elements from A and B into memory
				float64x2x4_t ABlocks = vld4q_f64(ARow + j * BlockSize);
				float64x2x4_t BBlocks = vld4q_f64(BRow + j * BlockSize);
				for (size_t ValIdx = 0; ValIdx < BlockSize / 2; ValIdx++)
				{
					ResBlocks.val[ValIdx] = vcmlaq_f64(ResBlocks.val[ValIdx], ABlocks.val[ValIdx], BBlocks.val[ValIdx]);
				}
				vst4q_f64(ResultRow + j * BlockSize, ResBlocks);
			}
			for (size_t j = B->NCols - Remainder; j < B->NCols; j++)
			{
				ResultRow[j] += ARow[j] * BRow[j];
			}
		}
	}
	return Result;
}

/*
v2 - Let's tile everything assuming we can completely pack our NEON registers.
We can do 4 FMLA/cycle on 128b-wide registers. So 8 float64s at a time.
The loops need to be swapped around a bit - we need to load the first block of A from outside the loop, and load multiple B blocks at once
Let's try 2x8 tiles

In debug mode this does indeed seem faster than the previous SIMD implementation... Still not as fast as the regular Savine version though.

Benchmarks on 1000x1000 matrix (these have some variability on each run, but the ranking stays the same):
Slow done in 3.470338s
Fast done in 1.222842s
FastER done in 2.462918s
FastEST done in 1.810831s

In release mode, it is faster!

Slow done in 0.170662s
Fast done in 0.156581s
FastER done in 0.357042s
FastEST done in 0.127912s

So maybe we just need to change the tile structure a bit. I don't think the float64x2x4_t structs actually block all the registers
we would like them to, it's not like they can execute all in parallel.
*/

void MultiplyTileSIMD2x8(Matrixf64 *Result, const Matrixf64 *A, const Matrixf64 *B, size_t TileX, size_t TileY)
{
	size_t TileRows = 2;
	size_t TileCols = 8;
	float64x2x4_t ResultRowBlocks[TileRows];
	float64x2x4_t ABlocks[TileRows];
	size_t ResultColIdx = TileY * TileCols;
	for (size_t i = 0; i < TileRows; i++)
	{
		size_t ResultRowIdx = TileX * TileRows + i;
		ResultRowBlocks[i] = vld4q_f64(Result->Data + ResultRowIdx * Result->NCols + ResultColIdx);
		ABlocks[i] = vld4q_f64(A->Data + ResultRowIdx * A->NCols + ResultColIdx);
	}
	for (size_t j = 0; j < Result->NRows; j++)
	{
		float64x2x4_t BBlocks = vld4q_f64(B->Data + j * B->NCols + ResultColIdx);
		for (size_t ValIdx = 0; ValIdx < 4; ValIdx++)
		{
			for (size_t i = 0; i < TileRows; i++)
			{
				ResultRowBlocks[i].val[ValIdx] = vcmlaq_f64(ResultRowBlocks[i].val[ValIdx], ABlocks[i].val[ValIdx], BBlocks.val[ValIdx]);
			}
		}
	}
	for (size_t i = 0; i < TileRows; i++)
	{
		size_t ResultRowIdx = TileX * TileRows + i;
		vst4q_f64(Result->Data + ResultRowIdx * Result->NCols + ResultColIdx, ResultRowBlocks[i]);
	}
}

/*
v3 - Let's do this in 2x4 mode now. Maybe a more efficient tiling structure

We could probably write a more generic implementation of this, but it's late and I'm lazy.

Verdict: In debug mode, this does not seem to make that much of a difference.... And it seems slower in release mode.

Benchmarks on 1000x1000:

Debug
Slow done in 3.451018s
Fast done in 1.214745s
FastER done in 2.437801s
FastEST 2x8 done in 1.821626s
FastEST 2x4 done in 1.880704s


Release
Slow done in 0.173660s
Fast done in 0.160243s
FastER done in 0.364735s
FastEST 2x8 done in 0.136150s
FastEST 2x4 done in 0.179667s
*/
void MultiplyTileSIMD2x4(Matrixf64 *Result, const Matrixf64 *A, const Matrixf64 *B, size_t TileX, size_t TileY)
{
	size_t TileRows = 2;
	size_t TileCols = 4;
	float64x2x2_t ResultRowBlocks[TileRows];
	float64x2x2_t ABlocks[TileRows];
	size_t ResultColIdx = TileY * TileCols;
	for (size_t i = 0; i < TileRows; i++)
	{
		size_t ResultRowIdx = TileX * TileRows + i;
		ResultRowBlocks[i] = vld2q_f64(Result->Data + ResultRowIdx * Result->NCols + ResultColIdx);
		ABlocks[i] = vld2q_f64(A->Data + ResultRowIdx * A->NCols + ResultColIdx);
	}
	for (size_t j = 0; j < Result->NRows; j++)
	{
		float64x2x2_t BBlocks = vld2q_f64(B->Data + j * B->NCols + ResultColIdx);
		for (size_t ValIdx = 0; ValIdx < 2; ValIdx++)
		{
			for (size_t i = 0; i < TileRows; i++)
			{
				ResultRowBlocks[i].val[ValIdx] = vcmlaq_f64(ResultRowBlocks[i].val[ValIdx], ABlocks[i].val[ValIdx], BBlocks.val[ValIdx]);
			}
		}
	}
	for (size_t i = 0; i < TileRows; i++)
	{
		size_t ResultRowIdx = TileX * TileRows + i;
		vst2q_f64(Result->Data + ResultRowIdx * Result->NCols + ResultColIdx, ResultRowBlocks[i]);
	}
}

void MultiplyTileSIMD1x8(Matrixf64 *Result, const Matrixf64 *A, const Matrixf64 *B, size_t TileX, size_t TileY)
{
	size_t TileCols = 8;
	size_t ResultColIdx = TileY * TileCols;
	float64x2x4_t ResultRowBlock = vld4q_f64(Result->Data + TileX * Result->NCols + ResultColIdx);
	float64x2x4_t ABlocks = vld4q_f64(A->Data + TileX * A->NCols + ResultColIdx);
	for (size_t j = 0; j < Result->NRows; j++)
	{
		float64x2x4_t BBlocks = vld4q_f64(B->Data + j * B->NCols + ResultColIdx);
		for (size_t ValIdx = 0; ValIdx < 4; ValIdx++)
		{
			ResultRowBlock.val[ValIdx] = vcmlaq_f64(ResultRowBlock.val[ValIdx], ABlocks.val[ValIdx], BBlocks.val[ValIdx]);
		}
	}
	vst4q_f64(Result->Data + TileX * Result->NCols + ResultColIdx, ResultRowBlock);
}

void MultiplyTileSIMD1x4(Matrixf64 *Result, const Matrixf64 *A, const Matrixf64 *B, size_t TileX, size_t TileY)
{
	size_t TileCols = 4;
	size_t ResultColIdx = TileY * TileCols;
	float64x2x2_t ResultRowBlock = vld2q_f64(Result->Data + TileX * Result->NCols + ResultColIdx);
	float64x2x2_t ABlocks = vld2q_f64(A->Data + TileX * A->NCols + ResultColIdx);
	for (size_t j = 0; j < Result->NRows; j++)
	{
		float64x2x4_t BBlocks = vld4q_f64(B->Data + j * B->NCols + ResultColIdx);
		for (size_t ValIdx = 0; ValIdx < 2; ValIdx++)
		{
			ResultRowBlock.val[ValIdx] = vcmlaq_f64(ResultRowBlock.val[ValIdx], ABlocks.val[ValIdx], BBlocks.val[ValIdx]);
		}
	}
	vst2q_f64(Result->Data + TileX * Result->NCols + ResultColIdx, ResultRowBlock);
}

/*
v4 Ok last try - let's load everything as 1x16 -> Maybe we get faster load times

Verdict.. Faster in debug mode. Not quite there in release mode.

1000x1000 Benchmarks

Debug
Slow done in 3.520341s
Fast done in 1.254194s
FastER done in 2.502608s
FastEST 2x8 done in 1.885986s
FastEST 2x4 done in 1.952919s
FastEST 1x16 done in 1.678129s

Slow done in 0.176381s
Fast done in 0.158877s
FastER done in 0.359372s
FastEST 2x8 done in 0.144498s
FastEST 2x4 done in 0.173988s
FastEST 1x16 done in 0.171464s
*/

void MultiplyTileSIMD1x16(Matrixf64 *Result, const Matrixf64 *A, const Matrixf64 *B, size_t TileX, size_t TileY)
{
	size_t TileCols = 16;
	size_t ResultColIdx = TileY * TileCols;
	float64x2x4_t ResultRowBlocks[2];
	float64x2x4_t ABlocks[2];
	ResultRowBlocks[0] = vld4q_f64(Result->Data + TileX * Result->NCols + ResultColIdx);
	ResultRowBlocks[1] = vld4q_f64(Result->Data + TileX * Result->NCols + ResultColIdx + 8);
	ABlocks[0] = vld4q_f64(A->Data + TileX * A->NCols + ResultColIdx);
	ABlocks[1] = vld4q_f64(A->Data + TileX * A->NCols + ResultColIdx + 8);
	float64x2x4_t BBlocks[2];
	for (size_t j = 0; j < Result->NRows; j++)
	{
		BBlocks[0] = vld4q_f64(B->Data + j * B->NCols + ResultColIdx);
		BBlocks[1] = vld4q_f64(B->Data + j * B->NCols + ResultColIdx + 8);
		for (size_t ValIdx = 0; ValIdx < 4; ValIdx++)
		{
			ResultRowBlocks[0].val[ValIdx] = vcmlaq_f64(ResultRowBlocks[0].val[ValIdx], ABlocks[0].val[ValIdx], BBlocks[0].val[ValIdx]);
			ResultRowBlocks[1].val[ValIdx] = vcmlaq_f64(ResultRowBlocks[1].val[ValIdx], ABlocks[1].val[ValIdx], BBlocks[1].val[ValIdx]);
		}
	}
	vst4q_f64(Result->Data + TileX * Result->NCols + ResultColIdx, ResultRowBlocks[0]);
	vst4q_f64(Result->Data + TileX * Result->NCols + ResultColIdx + 8, ResultRowBlocks[1]);
}

/*
v5 Of course! This should really be symetrical, so the "natural" tiling is 4x4. Let's give that a try

Verdict: In debug mode this does seem much better. Memory layout apparently is not quite optimal, as "Fast" still marginally beats us.

Release mode is fucking epic though, x2 speedup over Fast!

1000x1000 Benchmarks:

Debug:
Slow done in 3.495674s
Fast done in 1.255164s
FastER done in 2.502835s
FastEST 2x8 done in 1.890493s
FastEST 2x4 done in 1.918960s
FastEST 1x16 done in 3.247510s
FastEST 4x4 done in 1.512352s

Release:
Slow done in 0.171948s
Fast done in 0.153316s
FastER done in 0.356767s
FastEST 2x8 done in 0.099626s
FastEST 2x4 done in 0.163341s
FastEST 1x16 done in 0.238220s
FastEST 4x4 done in 0.086331s

Compare this with a single threaded numpy (OpenBLAS backend) benchmark of
0.042754173278808594s

We're starting to look decent
*/

void MultiplyTileSIMD4x4(Matrixf64 *Result, const Matrixf64 *A, const Matrixf64 *B, size_t TileX, size_t TileY)
{
	size_t TileRows = 4;
	size_t TileCols = 4;
	float64x2x2_t ResultRowBlocks[TileRows];
	float64x2x2_t ABlocks[TileRows];
	size_t ResultColIdx = TileY * TileCols;
	for (size_t i = 0; i < TileRows; i++)
	{
		size_t ResultRowIdx = TileX * TileRows + i;
		ResultRowBlocks[i] = vld2q_f64(Result->Data + ResultRowIdx * Result->NCols + ResultColIdx);
		ABlocks[i] = vld2q_f64(A->Data + ResultRowIdx * A->NCols + ResultColIdx);
	}
	for (size_t j = 0; j < Result->NRows; j++)
	{
		float64x2x2_t BBlocks = vld2q_f64(B->Data + j * B->NCols + ResultColIdx);
		for (size_t ValIdx = 0; ValIdx < 2; ValIdx++)
		{
			for (size_t i = 0; i < TileRows; i++)
			{
				ResultRowBlocks[i].val[ValIdx] = vcmlaq_f64(ResultRowBlocks[i].val[ValIdx], ABlocks[i].val[ValIdx], BBlocks.val[ValIdx]);
			}
		}
	}
	for (size_t i = 0; i < TileRows; i++)
	{
		size_t ResultRowIdx = TileX * TileRows + i;
		vst2q_f64(Result->Data + ResultRowIdx * Result->NCols + ResultColIdx, ResultRowBlocks[i]);
	}
}

void MultiplyTileSIMD4x8(Matrixf64 *Result, const Matrixf64 *A, const Matrixf64 *B, size_t TileX, size_t TileY)
{
	size_t TileRows = 4;
	size_t TileCols = 8;
	float64x2x4_t ResultRowBlocks[TileRows];
	float64x2x4_t ABlocks[TileRows];
	size_t ResultColIdx = TileY * TileCols;
	for (size_t i = 0; i < TileRows; i++)
	{
		size_t ResultRowIdx = TileX * TileRows + i;
		ResultRowBlocks[i] = vld4q_f64(Result->Data + (ResultRowIdx)*Result->NCols + ResultColIdx);
		ABlocks[i] = vld4q_f64(A->Data + (ResultRowIdx)*A->NCols + ResultColIdx);
		for (size_t j = 0; j < Result->NRows; j++)
		{
			float64x2x4_t BBlocks = vld4q_f64(B->Data + j * B->NCols + ResultColIdx);
			for (size_t ValIdx = 0; ValIdx < 4; ValIdx++)
			{
				ResultRowBlocks[i].val[ValIdx] = vcmlaq_f64(ResultRowBlocks[i].val[ValIdx], ABlocks[i].val[ValIdx], BBlocks.val[ValIdx]);
			}
		}
		vst4q_f64(Result->Data + ResultRowIdx * Result->NCols + ResultColIdx, ResultRowBlocks[i]);
	}
}

/*
v6
Now let's actually try using register blocking

1000x100 Benchmark - single run

Debug
Slow done in 3.448017s
Fast done in 1.213292s
FastER done in 2.427642s
FastEST 2x8 done in 1.874457s
FastEST 2x4 done in 1.905910s
FastEST 1x16 done in 1.644576s
FastEST 4x4 done in 1.543720s
FastEST 4x8 done in 1.900364s
MultiplyRegisterBlock 2x4 done in 2.300488s
MultiplyRegisterBlock 4x2 done in 2.282519s
MultiplyRegisterBlock 4x4 done in 2.258337s

Release
Slow done in 0.176373s
Fast done in 0.157340s
FastER done in 0.358622s
FastEST 2x8 done in 0.143927s
FastEST 2x4 done in 0.172985s
FastEST 1x16 done in 0.165954s
FastEST 4x4 done in 0.107556s
FastEST 4x8 done in 0.182878s
MultiplyRegisterBlock 2x4 done in 0.275106s
MultiplyRegisterBlock 4x2 done in 0.318297s
MultiplyRegisterBlock 4x4 done in 0.250176s

It's better than the naive NEON implementation... But I probably need to optimize the
loading pattern for the registers used, which means loading everything in order from A to B
rather than thinking of it as computing one patch of C at a time.
*/

void MultiplyRegisterBlocks(Matrixf64 *Result, Matrixf64 *A, Matrixf64 *B, size_t NRegA, size_t NRegB, size_t ARowIdx, size_t BColIdx)
{
	float64x2x4_t ResultBlocks[NRegA];
	float64x2x4_t BBlock;
	float64x2x4_t ABroadcast;
	for (size_t i = 0; i < NRegA; i++)
	{
		ResultBlocks[i] = vld4q_f64(Result->Data + (ARowIdx + i) * Result->NCols + BColIdx);
	}
	for (size_t ARegIdx = 0; ARegIdx < NRegA; ARegIdx++)
	{
		for (size_t BRowIdx = 0; BRowIdx < B->NRows; BRowIdx += NRegB)
		{
			size_t LastRegister = (BRowIdx + NRegB < B->NRows) ? NRegB : B->NRows - BRowIdx;
			for (size_t i = 0; i < LastRegister; i++)
			{
				size_t AColIdx = BRowIdx + i;
				BBlock = vld4q_f64(B->Data + (BRowIdx + i) * B->NCols + BColIdx);
				ABroadcast = vld4q_dup_f64(A->Data + (ARowIdx + ARegIdx) * A->NCols + AColIdx);
				for (size_t j = 0; j < 4; j++)
				{
					ResultBlocks[ARegIdx].val[j] = vcmlaq_f64(ResultBlocks[ARegIdx].val[j], ABroadcast.val[j], BBlock.val[j]);
				}
			}
		}
	}
	for (size_t i = 0; i < NRegA; i++)
	{
		vst4q_f64(Result->Data + (ARowIdx + i) * Result->NCols + BColIdx, ResultBlocks[i]);
	}
}

Matrixf64 Matrixf64MultiplyRegisterBlock(Matrixf64 *A, Matrixf64 *B, size_t NRegA, size_t NRegB)
{
	size_t RegSize = 8;
	assert(A->NCols == B->NRows);
	Matrixf64 Result;
	Result.NRows = A->NRows;
	Result.NCols = B->NCols;
	Result.Data = calloc(Result.NRows * Result.NCols, sizeof(double));
	size_t NRegX = A->NRows / NRegA;
	size_t LeftoverRegX = A->NRows % NRegA;
	size_t NRegY = B->NCols / RegSize;
	size_t LeftoverRegY = B->NCols % RegSize;
	for (size_t RegRowIdx = 0; RegRowIdx < NRegX; RegRowIdx++)
	{
		size_t ARow = RegRowIdx * NRegA;
		for (size_t RegColIdx = 0; RegColIdx < NRegY; RegColIdx++)
		{
			size_t BCol = RegColIdx * RegSize;
			MultiplyRegisterBlocks(&Result, A, B, NRegA, NRegB, ARow, BCol);
		}
	}

	for (size_t RowIdx = A->NRows - LeftoverRegX; RowIdx < A->NRows; RowIdx++)
	{
		for (size_t ColIdx = 0; ColIdx < B->NCols - LeftoverRegY; ColIdx += 8)
		{
			MultiplyTileSIMD1x8(&Result, A, B, RowIdx, ColIdx);
		}
	}
	for (size_t RowIdx = 0; RowIdx < A->NRows; RowIdx++)
	{
		const double *ARow = A->Data + RowIdx * A->NCols;
		double *ResultRow = Result.Data + RowIdx * Result.NCols;
		for (size_t k = 0; k < B->NRows; k++)
		{
			const double *BRow = B->Data + k * B->NCols;
			const double Aik = ARow[k];
			for (size_t ColIdx = B->NCols - LeftoverRegY; ColIdx < B->NCols; ColIdx++)
			{
				ResultRow[ColIdx] += Aik * BRow[ColIdx];
			}
		}
	}
	return Result;
}

/*
v7 - Let's actually load in order from A and B and accumulate the results back into C

1000x1000 Benchmarks

Debug
Slow done in 3.443199s
Fast done in 1.212221s
FastER done in 2.425718s
FastEST 2x8 done in 1.852756s
FastEST 2x4 done in 1.908038s
FastEST 1x16 done in 1.655259s
FastEST 4x4 done in 1.525683s
FastEST 4x8 done in 1.875789s
MultiplyRegisterBlock 2x4 done in 2.279805s
MultiplyRegisterBlock 4x2 done in 2.270494s
MultiplyRegisterBlock 4x4 done in 2.221423s
MultiplyAccumulate 4x4 done in 2.026345s
MultiplyAccumulate 2x4 done in 2.276662s
MultiplyAccumulate 4x2 done in 2.222507s

Release
Slow done in 0.174268s
Fast done in 0.155762s
FastER done in 0.356256s
FastEST 2x8 done in 0.129766s
FastEST 2x4 done in 0.169324s
FastEST 1x16 done in 0.162491s
FastEST 4x4 done in 0.100891s
FastEST 4x8 done in 0.180197s
MultiplyRegisterBlock 2x4 done in 0.268436s
MultiplyRegisterBlock 4x2 done in 0.313879s
MultiplyRegisterBlock 4x4 done in 0.254991s
MultiplyAccumulate 4x4 done in 0.305603s
MultiplyAccumulate 2x4 done in 0.397644s
MultiplyAccumulate 4x2 done in 0.394385s

Nope, this is even worse. Regular SIMD before the broadcasting was pretty good actually.
*/

void AccumulateRegisterBlocks(Matrixf64 *Result, double **AReg, double **BReg, size_t NRegA, size_t NRegB, size_t ARowIdx, size_t AColIdx, size_t BRowIdx, size_t BColIdx)
{
	float64x2x4_t ResultBlocks[NRegA];
	float64x2x4_t BBlock;
	float64x2x4_t ABroadcast;
	for (size_t i = 0; i < NRegA; i++)
	{
		ResultBlocks[i] = vld4q_f64(Result->Data + (ARowIdx + i) * Result->NCols + BColIdx);
	}
	for (size_t j = 0; j < NRegB; j++)
	{
		BBlock = vld4q_f64(BReg[j]);
		for (size_t i = 0; i < NRegA; i++)
		{
			ABroadcast = vld4q_dup_f64(AReg[i] + BRowIdx - AColIdx + j);
			for (size_t ValIdx = 0; ValIdx < 4; ValIdx++)
			{
				ResultBlocks[i].val[ValIdx] = vcmlaq_f64(ResultBlocks[i].val[ValIdx], ABroadcast.val[ValIdx], BBlock.val[ValIdx]);
			}
		}
	}
	for (size_t i = 0; i < NRegA; i++)
	{
		vst4q_f64(Result->Data + (ARowIdx + i) * Result->NCols + BColIdx, ResultBlocks[i]);
	}
}

// There is a corner case to figure out here (literally)
// Bottom left corner of a 1000x1000 matrix is two zeros with NRegA=NRegB=4... Why?
Matrixf64 Matrixf64MultiplyAccumulate(Matrixf64 *A, Matrixf64 *B, size_t NRegA, size_t NRegB)
{
	size_t RegSize = 8;
	assert(A->NCols == B->NRows);
	Matrixf64 Result;
	Result.NRows = A->NRows;
	Result.NCols = B->NCols;
	Result.Data = calloc(Result.NRows * Result.NCols, sizeof(double));
	size_t NRegX = A->NRows / NRegA;
	size_t ANRegCols = A->NCols / RegSize;
	size_t BNRowBlocks = RegSize / NRegB;
	size_t LeftoverRegX = A->NRows % NRegA;
	size_t NRegY = B->NCols / RegSize;
	size_t LeftoverRegY = B->NCols % RegSize;
	double *AReg[NRegA];
	double *BReg[NRegB];

	for (size_t ARegYIdx = 0; ARegYIdx < ANRegCols; ARegYIdx++)
	{
		size_t AColIdx = ARegYIdx * RegSize;
		for (size_t ARegBlock = 0; ARegBlock < NRegX; ARegBlock++)
		{
			size_t ARowIdx = ARegBlock * NRegA;
			for (size_t i = 0; i < NRegA; i++)
			{
				AReg[i] = A->Data + (ARowIdx + i) * A->NCols + AColIdx;
			}
			for (size_t BBlockRowIdx = 0; BBlockRowIdx < BNRowBlocks; BBlockRowIdx++)
			{
				size_t BRowIdx = AColIdx + BBlockRowIdx * NRegB;
				for (size_t BRegBlock = 0; BRegBlock < NRegY; BRegBlock++)
				{
					size_t BColIdx = BRegBlock * RegSize;
					for (size_t i = 0; i < NRegB; i++)
					{
						BReg[i] = B->Data + (BRowIdx + i) * B->NCols + BColIdx;
					}
					AccumulateRegisterBlocks(&Result, AReg, BReg, NRegA, NRegB, ARowIdx, AColIdx, BRowIdx, BColIdx);
				}
			}
		}
	}
	for (size_t RowIdx = A->NRows - LeftoverRegX; RowIdx < A->NRows; RowIdx++)
	{
		for (size_t ColIdx = 0; ColIdx < B->NCols - LeftoverRegY; ColIdx += 8)
		{
			MultiplyTileSIMD1x8(&Result, A, B, RowIdx, ColIdx);
		}
	}
	for (size_t RowIdx = 0; RowIdx < A->NRows; RowIdx++)
	{
		const double *ARow = A->Data + RowIdx * A->NCols;
		double *ResultRow = Result.Data + RowIdx * Result.NCols;
		for (size_t k = 0; k < B->NRows; k++)
		{
			const double *BRow = B->Data + k * B->NCols;
			const double Aik = ARow[k];
			for (size_t ColIdx = B->NCols - LeftoverRegY; ColIdx < B->NCols; ColIdx++)
			{
				ResultRow[ColIdx] += Aik * BRow[ColIdx];
			}
		}
	}
	return Result;
}

typedef enum tile_mode
{
	_2x8,
	_2x4,
	_1x16,
	_4x4,
	_4x8
} tile_mode;

Matrixf64 Matrixf64MultiplyFastest(const Matrixf64 *A, const Matrixf64 *B, tile_mode TileMode)
{
	assert(A->NCols == B->NRows);
	Matrixf64 Result;
	Result.NRows = A->NRows;
	Result.NCols = B->NCols;
	Result.Data = calloc(Result.NRows * Result.NCols, sizeof(double));
	size_t TileRows;
	size_t TileCols;
	switch (TileMode)
	{
	case _2x4:
		TileRows = 2;
		TileCols = 4;
		break;
	case _2x8:
		TileRows = 2;
		TileCols = 8;
		break;
	case _1x16:
		TileRows = 1;
		TileCols = 16;
	case _4x4:
		TileRows = 4;
		TileCols = 4;
		break;
	case _4x8:
		TileRows = 4;
		TileCols = 8;
		break;
	}
	size_t NTileRows = A->NRows / TileRows;
	size_t NLeftoverRows = A->NRows % TileRows;
	size_t NTileCols = B->NCols / TileCols;
	size_t NLeftoverCols = B->NCols % TileCols;
	for (size_t TileX = 0; TileX < NTileRows; TileX++)
	{
		for (size_t TileY = 0; TileY < NTileCols; TileY++)
		{
			switch (TileMode)
			{
			case _2x4:
				MultiplyTileSIMD2x4(&Result, A, B, TileX, TileY);
				break;
			case _2x8:
				MultiplyTileSIMD2x8(&Result, A, B, TileX, TileY);
				break;
			case _1x16:
				MultiplyTileSIMD1x16(&Result, A, B, TileX, TileY);
				break;
			case _4x4:
				MultiplyTileSIMD4x4(&Result, A, B, TileX, TileY);
				break;
			case _4x8:
				MultiplyTileSIMD4x8(&Result, A, B, TileX, TileY);
			}
		}
	}
	for (size_t RowIdx = A->NRows - NLeftoverRows; RowIdx < A->NRows; RowIdx++)
	{
		for (size_t TileY = 0; TileY < NTileCols; TileY++)
		{
			switch (TileMode)
			{
			case _2x4:
				MultiplyTileSIMD1x4(&Result, A, B, RowIdx, TileY);
				break;
			case _2x8:
				MultiplyTileSIMD1x8(&Result, A, B, RowIdx, TileY);
				break;
			case _1x16:
				break; // we've already done all the rows at this point
			case _4x4:
				MultiplyTileSIMD1x4(&Result, A, B, RowIdx, TileY);
				break;
			case _4x8:
				MultiplyTileSIMD1x8(&Result, A, B, RowIdx, TileY);
				break;
			}
		}
	}
	// Let's forgo SIMD on the last couple of columns because I'm lazy
	for (size_t RowIdx = 0; RowIdx < A->NRows; RowIdx++)
	{
		const double *ARow = A->Data + RowIdx * A->NCols;
		double *ResultRow = Result.Data + RowIdx * Result.NCols;
		for (size_t k = 0; k < B->NRows; k++)
		{
			const double *BRow = B->Data + k * B->NCols;
			const double Aik = ARow[k];
			for (size_t ColIdx = B->NCols - NLeftoverCols; ColIdx < B->NCols; ColIdx++)
			{
				ResultRow[ColIdx] += Aik * BRow[ColIdx];
			}
		}
	}
	return Result;
}

/* Naive implementation (that gets very optimized by the compiler) */

Matrixf64 Matrixf64MultiplyNaive(const Matrixf64 *A, const Matrixf64 *B)
{
	assert(A->NCols == B->NRows);
	Matrixf64 Result;
	Result.NRows = A->NRows;
	Result.NCols = B->NCols;
	Result.Data = calloc(Result.NRows * Result.NCols, sizeof(double));
	for (size_t i = 0; i < A->NRows; i++)
	{
		for (size_t j = 0; j < B->NCols; j++)
		{
			for (size_t k = 0; k < B->NRows; k++)
			{
				Result.Data[i * A->NRows + j] += A->Data[i * A->NRows + k] * B->Data[j * B->NRows + k];
			}
		}
	}
	return Result;
}

/*
Naive implementation of matrix multiplication. The thing is, there does not seem to be a noticeable difference
between this and the fast one above when turning on the -Ofast flag in clang.

The difference is noticeable at lower levels of optimization though. Some quick benchmarks for multiplying 2 random 1000 x 1000
matrices:

-O0:
Slow done in 3.454435s
Fast done in 1.221509s

-O1:
Slow done in 3.452416s
Fast done in 0.326448s

-O2:
Slow done in 0.841188s
Fast done in 0.158153s

-O3:
Slow done in 0.843141s
Fast done in 0.159574s

-Ofast:
Slow done in 0.163217s
Fast done in 0.156944s

Moral of the story, compilers are good at optimizing shit.

Numpy benchmark: running python 3.10, numpy 1.22.3, multiplying 2 1000x1000 random matrices is done in 0.0199429988861084s....

So python is still better. We need multithreading! With compiler optimizations on, we are already better than the numbers Savine reports
in his book as benchmarks on his Intel (on a MacBook Pro)

Pytorch 1.11.0:
Done in 0.0027608871459960938s

An order of magnitude faster than numpy?

*/

/*
Let's take the fast version and make it multithreaded

*/

typedef struct matrix_multiply_job
{
	double *ARow;
	double *BRow;
	size_t ARowIdx;
	size_t BRowIdx;
	size_t ARowLen;
	size_t BRowLen;
	volatile double *ResultRow;
	pthread_mutex_t *ResultRowLock;
} matrix_multiply_job;

typedef struct thread_queue
{
	size_t MaxQueueSize;
	size_t volatile NextEntryToDo;
	size_t volatile EntryCompletionCount;
	size_t volatile EntryCount;
	sem_t *Semaphore;
	pthread_mutex_t Lock;
} thread_queue;

void Matrixf64Multiply2Rows(matrix_multiply_job *Entry)
{
	pthread_mutex_lock(Entry->ResultRowLock);
	// #pragma clang loop vectorize(enable) // according to the compiler this fails even if we ask for it explicitely
	for (size_t k = 0; k < Entry->BRowLen; ++k)
	{
		// we just need to lock this so only one thread can write to it at a time
		Entry->ResultRow[k] += Entry->ARow[k] * Entry->BRow[k];
	}
	pthread_mutex_unlock(Entry->ResultRowLock);
}

void InitThreadQueue(thread_queue *Queue, size_t MaxQueueSize)
{
	pthread_mutex_init(&(Queue->Lock), NULL);
	Queue->Semaphore = sem_open("queue_semaphore", O_CREAT);
	Queue->MaxQueueSize = MaxQueueSize;
	Queue->NextEntryToDo = 0;
	Queue->EntryCount = 0;
	Queue->EntryCompletionCount = 0;
}

void PushContainerToQueue(thread_queue *Queue)
{
	Queue->EntryCount++;
	sem_post(Queue->Semaphore);
}

typedef struct thread_input
{
	thread_queue *Queue;
	matrix_multiply_job *Jobs;
	int ThreadIdx;
} thread_input;

bool IsWorkRemaining(thread_queue *Queue)
{
	pthread_mutex_lock(&(Queue->Lock));
	bool WorkLeft = Queue->EntryCompletionCount < Queue->EntryCount;
	pthread_mutex_unlock(&(Queue->Lock));
	return WorkLeft;
}

typedef struct work_queue_entry
{
	size_t EntryIdx;
	bool IsValid;
} work_queue_entry;

work_queue_entry GetNextElement(thread_queue *Queue)
{
	work_queue_entry Result;
	Result.IsValid = false;
	pthread_mutex_lock(&(Queue->Lock));
	if (Queue->NextEntryToDo < Queue->EntryCount)
	{
		Result.EntryIdx = Queue->NextEntryToDo++;
		Result.IsValid = true;
	}
	pthread_mutex_unlock(&(Queue->Lock));
	return Result;
}

void MarkEntryComplete(thread_queue *Queue)
{
	pthread_mutex_lock(&(Queue->Lock));
	Queue->EntryCompletionCount++;
	pthread_mutex_unlock(&(Queue->Lock));
}

bool DoMatrixMultiplyf64ThreadWork(matrix_multiply_job *Jobs, thread_queue *Queue)
{
	work_queue_entry Entry = GetNextElement(Queue);
	if (Entry.IsValid)
	{
		matrix_multiply_job *Job = Jobs + Entry.EntryIdx;
		Matrixf64Multiply2Rows(Job);
		MarkEntryComplete(Queue);
	}
	return Entry.IsValid;
}

void HandleMatrixMultiplyThread(void *Input)
{
	thread_input *InputArgs = (thread_input *)Input;
	for (;;)
	{
		if (IsWorkRemaining(InputArgs->Queue))
		{
			DoMatrixMultiplyf64ThreadWork(InputArgs->Jobs, InputArgs->Queue);
		}
		else
		{
			sem_wait(InputArgs->Queue->Semaphore);
		}
	}
}

Matrixf64
Matrixf64MultiplyMultithreaded(const Matrixf64 *A, const Matrixf64 *B, size_t NThreads)
{
	assert(A->NCols == B->NRows);
	thread_input ThreadInfo[NThreads];
	matrix_multiply_job *Jobs = malloc(sizeof(matrix_multiply_job) * A->NRows * B->NRows);
	pthread_t Threads[NThreads];
	thread_queue *Queue = malloc(sizeof(thread_queue));
	InitThreadQueue(Queue, 1000 * 1000 + 1);
	for (size_t i = 0; i < NThreads; i++)
	{
		ThreadInfo[i].ThreadIdx = i;
		ThreadInfo[i].Jobs = Jobs;
		ThreadInfo[i].Queue = Queue;
		pthread_create(Threads + i, NULL, (void *)&HandleMatrixMultiplyThread, (void *)(ThreadInfo + i));
	}
	Matrixf64 Result;
	Result.NRows = A->NRows;
	Result.NCols = B->NCols;
	Result.Data = calloc(Result.NRows * Result.NCols, sizeof(double));
	assert(A->NRows * B->NRows < Queue->MaxQueueSize);
	pthread_mutex_t RowLocks[A->NRows];
	for (size_t i = 0; i < A->NRows; i++)
	{
		pthread_mutex_init(RowLocks + i, NULL);
	}
	for (size_t j = 0; j < B->NRows; j++)
	{
		for (size_t i = 0; i < A->NRows; i++)
		{
			// go in this order so that we stack the jobs in column major of the result matrix -> maybe better order for the multithreading queue
			matrix_multiply_job *ThisJob = Jobs + j * B->NRows + i;
			ThisJob->ARowIdx = i;
			ThisJob->BRowIdx = j;
			ThisJob->ARow = A->Data + i * A->NCols;
			ThisJob->BRow = B->Data + j * B->NCols;
			ThisJob->ARowLen = A->NCols;
			ThisJob->BRowLen = B->NCols;
			ThisJob->ResultRow = Result.Data + i * A->NCols;
			ThisJob->ResultRowLock = RowLocks + i;
			PushContainerToQueue(Queue);
		}
	}
	while (IsWorkRemaining(Queue))
	{
		DoMatrixMultiplyf64ThreadWork(Jobs, Queue);
	}
	InitThreadQueue(Queue, 0); // set queue size back to 0 so all threads go back to sleep when funciton returns
	free(Jobs);
	return Result;
}

Matrixf64 RandomMatrixf64(size_t NRows, size_t NCols)
{
	srand(0);
	Matrixf64 Result;
	Result.NRows = NRows;
	Result.NCols = NCols;
	Result.Data = malloc(NRows * NCols * sizeof(double));
	for (size_t i = 0; i < NRows; i++)
	{
		for (size_t j = 0; j < NCols; j++)
		{
			Result.Data[i * NRows + j] = ((double)rand()) / RAND_MAX;
		}
	}
	return Result;
}

Matrixf64 Identity(size_t MatrixSize)
{
	Matrixf64 Result;
	Result.NRows = MatrixSize;
	Result.NCols = MatrixSize;
	Result.Data = calloc(MatrixSize * MatrixSize, sizeof(double));
	for (size_t i = 0; i < MatrixSize; i++)
	{
		Result.Data[i * MatrixSize + i] = 1;
	}
	return Result;
}

Matrixf64 Ones(size_t MatrixRows, size_t MatrixCols)
{
	Matrixf64 Result;
	Result.NRows = MatrixRows;
	Result.NCols = MatrixCols;
	Result.Data = malloc(MatrixRows * MatrixCols * sizeof(double));
	for (size_t i = 0; i < MatrixRows; i++)
	{
		for (size_t j = 0; j < MatrixCols; j++)
		{
			Result.Data[i * MatrixCols + j] = 1;
		}
	}
	return Result;
}

/* This is technically the slow version, but we get a x20 speedup with clang on -Ofast */
void test_matrix_multiply(size_t MatrixSize, int PrintResult, int NThreads, bool UseRandom)
{
	Matrixf64 A, B;
	if (UseRandom)
	{
		A = RandomMatrixf64(MatrixSize, MatrixSize);
		B = RandomMatrixf64(MatrixSize, MatrixSize);
	}
	else
	{
		A = Ones(MatrixSize, MatrixSize), B = Identity(MatrixSize);
	}
	double n_ops = MatrixSize * MatrixSize * MatrixSize * 2e-9;
	Time start = TimeGet();
	Matrixf64 Result = Matrixf64MultiplyNaive(&A, &B);
	Time end = TimeGet();
	free(Result.Data);
	printf("Slow done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	start = TimeGet();
	Result = Matrixf64MultiplyFast(&A, &B);
	end = TimeGet();
	printf("Fast done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	start = TimeGet();
	Result = Matrixf64MultiplyFaster(&A, &B);
	end = TimeGet();
	printf("FastER done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	tile_mode TileMode = _2x8;
	start = TimeGet();
	Result = Matrixf64MultiplyFastest(&A, &B, TileMode);
	end = TimeGet();
	printf("FastEST 2x8 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	TileMode = _2x4;
	start = TimeGet();
	Result = Matrixf64MultiplyFastest(&A, &B, TileMode);
	end = TimeGet();
	printf("FastEST 2x4 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	TileMode = _1x16;
	start = TimeGet();
	Result = Matrixf64MultiplyFastest(&A, &B, TileMode);
	end = TimeGet();
	printf("FastEST 1x16 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	TileMode = _4x4;
	start = TimeGet();
	Result = Matrixf64MultiplyFastest(&A, &B, TileMode);
	end = TimeGet();
	printf("FastEST 4x4 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	TileMode = _4x8;
	start = TimeGet();
	Result = Matrixf64MultiplyFastest(&A, &B, TileMode);
	end = TimeGet();
	printf("FastEST 4x8 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	size_t NRegA = 2, NRegB = 4;
	start = TimeGet();
	Result = Matrixf64MultiplyRegisterBlock(&A, &B, NRegA, NRegB);
	end = TimeGet();
	printf("MultiplyRegisterBlock 2x4 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	NRegA = 4;
	NRegB = 2;
	start = TimeGet();
	Result = Matrixf64MultiplyRegisterBlock(&A, &B, NRegA, NRegB);
	end = TimeGet();
	printf("MultiplyRegisterBlock 4x2 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	NRegA = 4;
	NRegB = 4;
	start = TimeGet();
	Result = Matrixf64MultiplyRegisterBlock(&A, &B, NRegA, NRegB);
	end = TimeGet();
	printf("MultiplyRegisterBlock 4x4 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	start = TimeGet();
	Result = Matrixf64MultiplyAccumulate(&A, &B, NRegA, NRegB);
	end = TimeGet();
	printf("MultiplyAccumulate 4x4 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	// start = TimeGet();
	free(Result.Data);
	NRegA = 2;
	NRegB = 4;
	start = TimeGet();
	Result = Matrixf64MultiplyAccumulate(&A, &B, NRegA, NRegB);
	end = TimeGet();
	printf("MultiplyAccumulate 2x4 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	free(Result.Data);
	NRegA = 4;
	NRegB = 2;
	start = TimeGet();
	Result = Matrixf64MultiplyAccumulate(&A, &B, NRegA, NRegB);
	end = TimeGet();
	printf("MultiplyAccumulate 4x2 done in %lf\n", n_ops / TimeDiffInSeconds(start, end));
	// free(Result.Data);
	// Result = Matrixf64MultiplyMultithreaded(&A, &B, NThreads);
	// end = TimeGet();
	// printf("Multithreaded done in %lfs with %d threads\n", n_ops / TimeDiffInSeconds(start, end), NThreads);
	if (PrintResult)
	{
		printf("Result\n");
		for (size_t i = 0; i < Result.NRows; i++)
		{
			for (size_t j = 0; j < Result.NCols; j++)
			{
				printf("%f ", Result.Data[i * Result.NRows + j]);
			}
			printf("\n");
		}
	}
}

int main()
{
	test_matrix_multiply(64, 0, 2, 0);
	return 0;
}
#endif
