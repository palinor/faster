#pragma once
#ifdef __APPLE__
#include <arm_neon.h>
#else 
#include <immintrin.h>
#endif

#include <math.h>
#include <cassert>
#include <cstdlib>

#include "threadpool.cpp"

#define MAX_FLOAT 1e20
#define MIN_FLOAT 1e-12

struct matrix_f32
{
	size_t rows;
	size_t cols;
	float *contents;
	size_t ld;
	float *operator[](int i) { return &contents[i * ld]; }
};


inline void FreeMatrixf32(matrix_f32 *a) {
	free(a->contents);
}

inline float MatrixGetItem(matrix_f32 *a, size_t i, size_t j)
{
	return a->contents[i * a->ld + j];
}

inline void MatrixPutItem(matrix_f32 *a, size_t i, size_t j, float x) {
	a->contents[i * a->ld + j] = x;
}

inline float *MatrixGetAddr(const matrix_f32 *a, size_t i, size_t j)
{
	return a->contents + (i * a->ld + j);
}

inline float FAbs(float x) {
	if (x < 0) {
		return -x;
	}
	else {
		return x;
	}
}

inline int IAbs(int x) {
	if (x < 0) {
		return -x;
	}
	else {
		return x;
	}
}

matrix_f32 Matrixf32Copy(matrix_f32 *a) {
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = a->cols;
	result.ld = a->ld;
	result.contents = (float *)malloc(a->rows * a->cols * sizeof(float));
	if (!result.contents) {
		perror("Error allocating result.contents in Matrixf32Copy");
		exit(1);
	}
	for (size_t i = 0; i < a->rows; i++) {
		for (size_t j = 0; j < a->cols; j++) {
			MatrixPutItem(&result, i, j, MatrixGetItem(a, i, j));
		}
	}
	return result;
}



void PrintMatrixf32(matrix_f32 *result)
{
	for (size_t i = 0; i < result->rows; i++)
	{
		for (size_t j = 0; j < result->cols; j++)
		{
			printf("%f ", MatrixGetItem(result, i, j));
		}
		printf("\n");
	}
}

double GetMaxError(matrix_f32 *a, matrix_f32 *b)
{
	assert(a->rows == b->rows);
	assert(a->cols == b->cols);
	float maxError = 0;
	for (size_t i = 0; i < a->rows; i++)
	{
		for (size_t j = 0; j < a->cols; j++)
		{
			float thisError = FAbs(MatrixGetItem(a, i, j) - MatrixGetItem(b, i, j));
			if (!(thisError < MAX_FLOAT)) {
				return MAX_FLOAT; // handle NaN
			}
			maxError = (thisError > maxError) ? thisError : maxError;
		}
	}
	return maxError;
}


int Matrixf32MultiplyToTarget(matrix_f32 *result, const matrix_f32 *a, const matrix_f32 *b)
{
	assert(a->cols == b->rows);
	// We check that result has already been allocated properly
	assert(result->rows == a->rows);
	assert(result->cols == b->cols);
	assert(result->ld == b->cols);
	assert(!!(result->contents));
	if (!result->contents) {
		perror("Error allocating result.contents in Matrixf32Multiply");
		exit(1);
	}
	for (size_t i = 0; i < a->rows; i++)
	{
		const float *aRow = a->contents + i * a->ld;
		float *resultRow = result->contents + i * result->ld;
		for (size_t k = 0; k < b->rows; k++)
		{
			const float *bRow = b->contents + k * b->ld;
			const float Aik = aRow[k];
			// #pragma clang loop vectorize_width(8) // interleave_count(4)
			for (size_t j = 0; j < result->cols; j++)
			{
				//todo(AION): Have this swap depending on the platform we compile on.
				//resultRow[j] += Aik * bRow[j];
				/*
				This part is actually out of line with the microkernel loop (no difference in -Ofast, but visible in other compile modes)
				The above lines up just fine in all

				BUUUUUUUT.... not on windows. So we actually need to swap depending on the compiler. Lol.
				*/
				resultRow[j] = fma(bRow[j], Aik, resultRow[j]);
			}
		}
	}
	return 0;
}

matrix_f32 Matrixf32Multiply(const matrix_f32 *a, const matrix_f32 *b)
{
	assert(a->cols == b->rows);
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = b->cols;
	result.contents = (float *)calloc(result.rows * result.cols, sizeof(float));
	Matrixf32MultiplyToTarget(&result, a, b);
	return result;
}

/*
Multiply inplace by a scalar
*/
void Matrixf32MultiplyInplace(matrix_f32 *a, float x) {
	for (size_t i = 0; i < a->rows; i++) {
		for (size_t j = 0; j < a->cols; j++) {
			MatrixPutItem(a, i, j, x * MatrixGetItem(a, i, j));
		}
	}
}

matrix_f32 Matrixf32Sum(matrix_f32 *a, matrix_f32 *b) {
	assert(a->rows == b->rows);
	assert(a->cols == b->cols);
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = result.cols;
	size_t contentSize = result.rows * result.cols * sizeof(float);
	void *newContent = malloc(contentSize);
	if (!newContent) {
		perror("Error allocating results in Matrixf32Sum");
		exit(1);
	}
	result.contents = (float *)newContent;
	for (size_t i = 0; i < result.rows; i++) {
		for (size_t j = 0; j < result.cols; j++) {
			MatrixPutItem(&result, i, j, MatrixGetItem(a, i, j) + MatrixGetItem(b, i, j));
		}
	}
	return result;
}


matrix_f32 RandomMatrixf32(size_t rows, size_t cols)
{
	matrix_f32 result;
	result.rows = rows;
	result.cols = cols;
	result.ld = cols;
	void *newContents = malloc(rows * cols * sizeof(float));
	if (!newContents) {
		perror("Error allocating result.contents in RandomMatrixf32");
		exit(1);
	}
	result.contents = (float *)newContents;
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			int randInt = rand();
			float randomFloat = ((float)randInt) / RAND_MAX;
			MatrixPutItem(&result, i, j, randomFloat);
		}
	}
	return result;
}

matrix_f32 Onesf32(size_t rows, size_t cols)
{
	matrix_f32 result;
	result.rows = rows;
	result.cols = cols;
	result.ld = cols;
	result.contents = (float *)malloc(rows * cols * sizeof(float));
	if (!result.contents) {
		perror("Error allocating result.contents in Onesf32");
		exit(1);
	}

	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			MatrixPutItem(&result, i, j, 1);
		}
	}
	return result;
}

matrix_f32 Identityf32(size_t n)
{
	matrix_f32 result;
	result.rows = n;
	result.cols = n;
	result.ld = n;
	result.contents = (float *)calloc(n * n, sizeof(float));
	for (size_t i = 0; i < n; i++)
	{
		MatrixPutItem(&result, i, i, 1);
	}
	return result;
}

/*
Returns an upper triangular matrix with random non-zero elements
*/
matrix_f32 UpperRandomf32(size_t n) {
	matrix_f32 result;
	result.rows = n;
	result.cols = n;
	result.ld = n;
	result.contents = (float *)calloc(n * n, sizeof(float));
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i; j < n; j++) {
			int randInt = rand();
			float randomFloat = ((float)randInt) / RAND_MAX;
			MatrixPutItem(&result, i, j, randomFloat * 100);
		}
	}
	return result;
}

/*
Returns an lower triangular matrix with random non-zero elements
*/
matrix_f32 LowerRandomf32(size_t n) {
	matrix_f32 result;
	result.rows = n;
	result.cols = n;
	result.ld = n;
	result.contents = (float *)calloc(n * n, sizeof(float));
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j > 0; j--) { // stupid trick because j is an unsigned int
			int randomValue = rand();
			float newValue = static_cast<float>(randomValue) / RAND_MAX * 100;
			MatrixPutItem(&result, i, j - 1, newValue);
		}
	}
	return result;
}

/*
Internal block for matrix multiplication - compute the result in blocks of isze kernelHeight x (4 * kernelWidth)
*/
#ifdef _M_ARM
// ARM implimentation
void Microkernel(const matrix_f32 *a, const matrix_f32 *b, matrix_f32 *result, size_t resultRow, size_t resultCol, size_t kernelWidth, size_t kernelHeight)
{
	assert(a->cols == b->rows);
	float32x4_t resultElems[kernelHeight][kernelWidth];
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			resultElems[i][j] = vdupq_n_f32(0);
		}
	}
	float32x4_t aElems[kernelHeight];
	float32x4_t bElems[4][kernelWidth];
	for (size_t k = 0; k < a->cols / 4; k++)
	{
		for (size_t i = 0; i < kernelHeight; i++)
		{
			aElems[i] = vld1q_f32(MatrixGetAddr(a, resultRow + i, 4 * k));
			for (size_t j = 0; j < kernelWidth; j++)
			{
				bElems[0][j] = vld1q_f32(MatrixGetAddr(b, 4 * k, resultCol + 4 * j));
				resultElems[i][j] = vmlaq_laneq_f32(resultElems[i][j], bElems[0][j], aElems[i], 0);

				bElems[1][j] = vld1q_f32(MatrixGetAddr(b, 4 * k + 1, resultCol + 4 * j));
				resultElems[i][j] = vmlaq_laneq_f32(resultElems[i][j], bElems[1][j], aElems[i], 1);

				bElems[2][j] = vld1q_f32(MatrixGetAddr(b, 4 * k + 2, resultCol + 4 * j));
				resultElems[i][j] = vmlaq_laneq_f32(resultElems[i][j], bElems[2][j], aElems[i], 2);

				bElems[3][j] = vld1q_f32(MatrixGetAddr(b, 4 * k + 3, resultCol + 4 * j));
				resultElems[i][j] = vmlaq_laneq_f32(resultElems[i][j], bElems[3][j], aElems[i], 3);
			}
		}
	}

	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			vst1q_f32(MatrixGetAddr(result, resultRow + i, resultCol + 4 * j), resultElems[i][j]);
		}
	}
}

#else
// AMD64 implementation
#define SIMD_VECTOR_SIZE 8 
void Microkernel(matrix_f32 *resultMatrix, const matrix_f32 *a, const matrix_f32 *b, size_t resultRow, size_t resultCol, const size_t kernelWidth, const size_t kernelHeight)
{
	assert(a->cols == b->rows);
	assert(kernelHeight * kernelWidth < 32);
	__m256 resultElems[32];
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			resultElems[i * kernelWidth + j] = _mm256_set1_ps(0);
		}
	}
	for (size_t k = 0; k < a->cols / SIMD_VECTOR_SIZE; k++)
	{
		for (size_t i = 0; i < kernelHeight; i++)
		{
			float *aElems_i = MatrixGetAddr(a, resultRow + i, SIMD_VECTOR_SIZE * k);
			for (size_t j = 0; j < kernelWidth; j++)
			{
				//todo(AION): we have 3 fma ports, why can't we use 3 accumulators? It seemed needlessly slow when I tried 2
				__m256 bElems_0j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i0 = _mm256_set1_ps(*aElems_i);
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i0, bElems_0j, resultElems[i * kernelWidth + j]);

				__m256 bElems_1j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 1, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i1 = _mm256_set1_ps(*(aElems_i + 1));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i1, bElems_1j, resultElems[i * kernelWidth + j]);

				__m256 bElems_2j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 2, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i2 = _mm256_set1_ps(*(aElems_i + 2));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i2, bElems_2j, resultElems[i * kernelWidth + j]);

				__m256 bElems_3j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 3, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i3 = _mm256_set1_ps(*(aElems_i + 3));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i3, bElems_3j, resultElems[i * kernelWidth + j]);

				__m256 bElems_4j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 4, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i4 = _mm256_set1_ps(*(aElems_i + 4));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i4, bElems_4j, resultElems[i * kernelWidth + j]);

				__m256 bElems_5j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 5, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i5 = _mm256_set1_ps(*(aElems_i + 5));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i5, bElems_5j, resultElems[i * kernelWidth + j]);

				__m256 bElems_6j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 6, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i6 = _mm256_set1_ps(*(aElems_i + 6));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i6, bElems_6j, resultElems[i * kernelWidth + j]);

				__m256 bElems_7j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 7, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i7 = _mm256_set1_ps(*(aElems_i + 7));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i7, bElems_7j, resultElems[i * kernelWidth + j]);
			}
		}
	}
	// If k is not divisible by SIMD_VECTOR_SIZE, handle the remainder
	size_t remainder = a->cols % SIMD_VECTOR_SIZE;
	if (remainder) {
		for (size_t i = 0; i < kernelHeight; i++)
		{
			float *aRemainder_i = MatrixGetAddr(a, resultRow + i, a->cols - remainder);
			for (size_t j = 0; j < kernelWidth; j++)
			{
				__m256 bRemainderElems_0j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i0 = _mm256_set1_ps(*aRemainder_i);
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i0, bRemainderElems_0j, resultElems[i * kernelWidth + j]);
				if (remainder > 1)
				{
					__m256 bRemainderElems_1j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 1, resultCol + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i1 = _mm256_set1_ps(*(aRemainder_i + 1));
					resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i1, bRemainderElems_1j, resultElems[i * kernelWidth + j]);
				}
				if (remainder > 2)
				{
					__m256 bRemainderElems_2j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 2, resultCol + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i2 = _mm256_set1_ps(*(aRemainder_i + 2));
					resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i2, bRemainderElems_2j, resultElems[i * kernelWidth + j]);
				}
				if (remainder > 3)
				{
					__m256 bRemainderElems_3j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 3, resultCol + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i3 = _mm256_set1_ps(*(aRemainder_i + 3));
					resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i3, bRemainderElems_3j, resultElems[i * kernelWidth + j]);
				}
				if (remainder > 4)
				{
					__m256 bRemainderElems_4j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 4, resultCol + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i4 = _mm256_set1_ps(*(aRemainder_i + 4));
					resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i4, bRemainderElems_4j, resultElems[i * kernelWidth + j]);
				}
				if (remainder > 5)
				{
					__m256 bRemainderElems_5j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 5, resultCol + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i5 = _mm256_set1_ps(*(aRemainder_i + 5));
					resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i5, bRemainderElems_5j, resultElems[i * kernelWidth + j]);
				}
				if (remainder > 6)
				{
					__m256 bRemainderElems_6j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 6, resultCol + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i6 = _mm256_set1_ps(*(aRemainder_i + 6));
					resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i6, bRemainderElems_6j, resultElems[i * kernelWidth + j]);
				}
			}
		}
	}
	float *resultLocation = MatrixGetAddr(resultMatrix, resultRow, resultCol);
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			_mm256_store_ps(resultLocation + i * b->cols + j * SIMD_VECTOR_SIZE, resultElems[i * kernelWidth + j]);
		}
	}

}

/*
If k is not divisible by SIMD_VECTOR_SIZE, handle the addition of the remaining elements
void MicrokernelRemainder(const matrix_f32 *a, const matrix_f32 *b, float *resultLocation, size_t resultRow, size_t resultCol, size_t kernelWidth, size_t kernelHeight)
{
	size_t remainder = a->cols % SIMD_VECTOR_SIZE;
	__m256 resultElems[SIMD_VECTOR_SIZE * 4];
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			resultElems[i * kernelWidth + j] = _mm256_load_ps(0);
		}
	}
	for (size_t i = 0; i < kernelHeight; i++)
	{
		float *aRemainder_i = MatrixGetAddr(a, resultRow + i, a->cols - remainder);
		for (size_t j = 0; j < kernelWidth; j++)
		{
			__m256 bRemainderElems_0j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder, resultCol + SIMD_VECTOR_SIZE * j));
			__m256 aRemainderElems_i0 = _mm256_set1_ps(*aRemainder_i);
			resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i0, bRemainderElems_0j, resultElems[i * kernelWidth + j]);
			if (remainder > 1)
			{
				__m256 bRemainderElems_1j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 1, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i1 = _mm256_set1_ps(*(aRemainder_i + 1));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i1, bRemainderElems_1j, resultElems[i * kernelWidth + j]);
			}
			if (remainder > 2)
			{
				__m256 bRemainderElems_2j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 2, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i2 = _mm256_set1_ps(*(aRemainder_i + 2));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i2, bRemainderElems_2j, resultElems[i * kernelWidth + j]);
			}
			if (remainder > 3)
			{
				__m256 bRemainderElems_3j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 3, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i3 = _mm256_set1_ps(*(aRemainder_i + 3));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i3, bRemainderElems_3j, resultElems[i * kernelWidth + j]);
			}
			if (remainder > 4)
			{
				__m256 bRemainderElems_4j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 4, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i4 = _mm256_set1_ps(*(aRemainder_i + 4));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i4, bRemainderElems_4j, resultElems[i * kernelWidth + j]);
			}
			if (remainder > 5)
			{
				__m256 bRemainderElems_5j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 5, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i5 = _mm256_set1_ps(*(aRemainder_i + 5));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i5, bRemainderElems_5j, resultElems[i * kernelWidth + j]);
			}
			if (remainder > 6)
			{
				__m256 bRemainderElems_6j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 6, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i6 = _mm256_set1_ps(*(aRemainder_i + 6));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i6, bRemainderElems_6j, resultElems[i * kernelWidth + j]);
			}
		}
	}
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			_mm256_store_ps(resultLocation, resultElems[i * kernelWidth + j]);
		}
	}
}*/
#endif

struct matrix_microkernel_row_task_info {
	const matrix_f32 *leftMatrix;
	const matrix_f32 *rightMatrix;
	matrix_f32 *resultMatrix;
	size_t resultRowNumber;
	size_t kernelWidth;
	size_t kernelHeight;
};


void MatrixMicrokernelRowTask(void *taskInfoInput) {
	matrix_microkernel_row_task_info *taskInfo = reinterpret_cast<matrix_microkernel_row_task_info *>(taskInfoInput);
	const matrix_f32 *a = taskInfo->leftMatrix;
	const matrix_f32 *b = taskInfo->rightMatrix;
	size_t row = taskInfo->resultRowNumber;
	size_t kernelWidth = taskInfo->kernelWidth;
	size_t kernelHeight = taskInfo->kernelHeight;
	matrix_f32 *resultMatrix = taskInfo->resultMatrix;
	size_t nCols = taskInfo->rightMatrix->cols / (kernelWidth * SIMD_VECTOR_SIZE);
	for (size_t col = 0; col < nCols; col++) {
		Microkernel(resultMatrix, a, b, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, kernelHeight);
	}
	free(taskInfoInput);
}



void Matrixf32MicrokernelMultiply(
	matrix_f32 *result,
	const matrix_f32 *a,
	const matrix_f32 *b,
	const size_t kernelWidth,
	const size_t kernelHeight,
	thread_pool *pool
)
{
	assert(a->cols == b->rows);
	assert(result->rows == a->rows);
	assert(result->cols == b->cols);
	assert(!!result->contents);
	const size_t nRows = a->rows / kernelHeight;
	const size_t rowRemainder = a->rows % kernelHeight;
	const size_t nCols = b->cols / (kernelWidth * SIMD_VECTOR_SIZE);
	const size_t colRemainder = a->cols % SIMD_VECTOR_SIZE;
	// task_handle *futures = (task_handle *)malloc(sizeof(task_handle) * nRows);
	// Do the main block
	for (int row = 0; row < nRows; row++) {
		// Each row defines a task. These are then run in parallel if there is a thread pool
		// otherwise they are run sequentially 
		/*
		auto rowTask = [=, &result]() {
			for (size_t col = 0; col < nCols; col++) {
				float *resultLocation = MatrixGetAddr(&result, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE);
				Microkernel(a, b, resultLocation, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, kernelHeight);
				if (colRemainder) {
					MicrokernelRemainder(a, b, &result, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, kernelHeight);
				}
			}
			return true;
		}; */
		matrix_microkernel_row_task_info *taskInfo = (matrix_microkernel_row_task_info *)malloc(sizeof(matrix_microkernel_row_task_info));
		taskInfo->kernelWidth = kernelWidth;
		taskInfo->kernelHeight = kernelHeight;
		taskInfo->leftMatrix = a;
		taskInfo->rightMatrix = b;
		taskInfo->resultRowNumber = row;
		taskInfo->resultMatrix = result;


		if (pool) {
			PushTaskToQueue(taskInfo, &MatrixMicrokernelRowTask, &(pool->queue_));
			// task_handle rowFuture = PushTaskToQueue(rowTask, &(pool->queue_));
			// new(futures + row) task_handle(std::move(rowFuture));
		}
		else {
			MatrixMicrokernelRowTask(&taskInfo);
		}
	}
	if (pool) {
		ThreadPoolActiveWait(pool);
		ResetConcurrentTaskQueue(&(pool->queue_));
	}
	// Clean up what is left (columns)
	if (nCols * kernelWidth * SIMD_VECTOR_SIZE < b->cols)
	{
		for (size_t i = 0; i < a->rows; i++)
		{
			const float *aRow = a->contents + i * a->cols;
			float *resultRow = result->contents + i * result->cols;
			for (size_t k = 0; k < b->rows; k++)
			{
				const float *bRow = b->contents + k * b->cols;
				const float Aik = aRow[k];
				for (size_t j = nCols * kernelWidth * SIMD_VECTOR_SIZE; j < b->cols; j++)
				{
					//todo(AION): same as above, this needs to be compiler dependant to be in line
					//resultRow[j] += Aik * bRow[j];
					resultRow[j] = fma(bRow[j], Aik, resultRow[j]);

				}
			}
		}
	}
}

