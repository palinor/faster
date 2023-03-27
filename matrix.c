#include <arm_neon.h>
#include <assert.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#define ABS(x) (((x) > 0) ? (x) : -(x))
#define MAX_FLOAT 1e20
#define MIN_FLOAT 1e-12
#define MATRIX_ITEM(a, i, j) ((a)->contents[(i) * (a)->ld + (j)])
#define SIGN(x) (((x) > 0) ? 1 : -1)

typedef struct matrix_f32
{
	size_t rows;
	size_t cols;
	float *contents;
	size_t ld;
} matrix_f32;

void Free(matrix_f32 *a) {
	free(a->contents);
}


matrix_f32 Copy(matrix_f32 *a) {
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = a->cols;
	result.ld = a->ld;
	result.contents = malloc(a->rows * a->cols * sizeof(float));
	for (size_t i = 0; i < a->rows; i++) {
		for (size_t j = 0; j < a->cols; j++) {
			MATRIX_ITEM(&result, i, j) = MATRIX_ITEM(a, i, j);
		}
	}
	return result;
}

float MatrixGetItem(matrix_f32 *a, size_t i, size_t j)
{
	return a->contents[i * a->ld + j];
}

void MatrixPutItem(matrix_f32 *a, size_t i, size_t j, float x) {
	a->contents[i * a->ld + j] = x;
}

float *MatrixGetAddr(const matrix_f32 *a, size_t i, size_t j)
{
	return a->contents + i * a->ld + j;
}


void PrintMatrix(matrix_f32 *result)
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
			float thisError = ABS(MatrixGetItem(a, i, j) - MatrixGetItem(b, i, j));
			if (!(thisError < MAX_FLOAT)) {
				return MAX_FLOAT; // handle NaN
			}
			maxError = (thisError > maxError) ? thisError : maxError;
		}
	}
	return maxError;
}

matrix_f32 Multiply(const matrix_f32 *a, const matrix_f32 *b)
{
	assert(a->cols == b->rows);
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = b->cols;
	result.contents = calloc(result.rows * result.cols, sizeof(double));
	for (size_t i = 0; i < a->rows; i++)
	{
		const float *aRow = a->contents + i * a->cols;
		float *resultRow = result.contents + i * result.cols;
		for (size_t k = 0; k < b->rows; k++)
		{
			const float *bRow = b->contents + k * b->cols;
			const float Aik = aRow[k];
			// #pragma clang loop vectorize_width(8) // interleave_count(4)
			for (size_t j = 0; j < b->cols; j++)
			{
				resultRow[j] += Aik * bRow[j];
				/*
				This part is actually out of line with the microkernel loop (no difference in -Ofast, but visible in other compile modes)
				The above lines up just fine in all
				*/
				// resultRow[j] = fma(bRow[j], Aik, resultRow[j]);
			}
		}
	}
	return result;
}

/*
Multiply inplace by a scalar
*/
void MultiplyInplace(matrix_f32 *a, float x) {
	for (size_t i = 0; i < a->rows; i++) {
		for (size_t j = 0; j < a->cols; j++) {
			MATRIX_ITEM(a, i, j) *= x;
		}
	}
}

matrix_f32 Sum(matrix_f32 *a, matrix_f32 *b) {
	assert(a->rows == b->rows);
	assert(a->cols == b->cols);
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = result.cols;
	result.contents = malloc(a->rows * a->cols * sizeof(float));
	for (size_t i = 0; i < result.rows; i++) {
		for (size_t j = 0; j < result.cols; j++) {
			result.contents[i * result.ld + j] = MATRIX_ITEM(a, i, j) + MATRIX_ITEM(b, i, j);
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
	result.contents = malloc(rows * cols * sizeof(float));
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			result.contents[i * cols + j] = 1000 * ((float)rand()) / RAND_MAX;
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
	result.contents = malloc(rows * cols * sizeof(float));
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			result.contents[i * cols + j] = 1;
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
	result.contents = calloc(n * n, sizeof(float));
	for (size_t i = 0; i < n; i++)
	{
		result.contents[i * result.ld + i] = 1;
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
	result.contents = calloc(n * n, sizeof(float));
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i; j < n; j++) {
			result.contents[i * result.ld + j] = 100 * ((float)rand()) / RAND_MAX;
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
	result.contents = calloc(n * n, sizeof(float));
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j > 0; j--) { // stupid trick because j is an unsigned int
			MATRIX_ITEM(&result, i, j - 1) = 100 * ((float)rand()) / RAND_MAX;
		}
	}
	return result;
}

/*
Internal block for matrix multiplication - compute the result in blocks of isze kernelHeight x (4 * kernelWidth)
*/
void Microkernel(const matrix_f32 *a, const matrix_f32 *b, matrix_f32 *result, size_t resultRow, size_t resultCol, size_t kernelHeight, size_t kernelWidth)
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

/*
If k is not divisible by 4, handle the addition of the remaining elements
*/
void MicrokernelRemainder(const matrix_f32 *a, const matrix_f32 *b, matrix_f32 *result, size_t resultRow, size_t resultCol, size_t kernelHeight, size_t kernelWidth)
{
	size_t remainder = a->cols % 4;
	float32x4_t resultElems[kernelHeight][kernelWidth];
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			resultElems[i][j] = vld1q_f32(MatrixGetAddr(result, resultRow + i, resultCol + 4 * j));
		}
	}
	float32x4_t aRemainderElems[kernelHeight];
	float32x4_t bRemainderElems[remainder][kernelWidth];
	for (size_t i = 0; i < kernelHeight; i++)
	{
		aRemainderElems[i] = vld1q_f32(MatrixGetAddr(a, resultRow + i, a->cols - remainder));
		vsetq_lane_f32(0, aRemainderElems[i], 3);
		if (remainder < 3) {
			vsetq_lane_f32(0, aRemainderElems[i], 2);
		}
		if (remainder < 2) {
			vsetq_lane_f32(0, aRemainderElems[i], 1);
		}
		for (size_t j = 0; j < kernelWidth; j++)
		{
			bRemainderElems[0][j] = vld1q_f32(MatrixGetAddr(b, a->cols - remainder, resultCol + 4 * j));
			resultElems[i][j] = vmlaq_laneq_f32(resultElems[i][j], bRemainderElems[0][j], aRemainderElems[i], 0);
			if (remainder > 1)
			{
				bRemainderElems[1][j] = vld1q_f32(MatrixGetAddr(b, a->cols - remainder + 1, resultCol + 4 * j));
				resultElems[i][j] = vmlaq_laneq_f32(resultElems[i][j], bRemainderElems[1][j], aRemainderElems[i], 1);
			}
			if (remainder > 2)
			{
				bRemainderElems[2][j] = vld1q_f32(MatrixGetAddr(b, a->cols - remainder + 2, resultCol + 4 * j));
				resultElems[i][j] = vmlaq_laneq_f32(resultElems[i][j], bRemainderElems[2][j], aRemainderElems[i], 2);
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

matrix_f32 MatrixMicrokernelMultiply(const matrix_f32 *a, const matrix_f32 *b)
{
	assert(a->cols == b->rows);
	size_t kernelHeight = 4;
	size_t kernelWidth = 4;
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = b->cols;
	result.contents = calloc(a->rows * b->cols, sizeof(float));
	size_t nRows = a->rows / kernelHeight;
	size_t nCols = b->cols / (kernelWidth * 4);
	size_t remainder = a->cols % 4;
	// Do the main block
	for (size_t row = 0; row < nRows; row++)
	{
		for (size_t col = 0; col < nCols; col++)
		{
			Microkernel(a, b, &result, row * kernelHeight, col * kernelWidth * 4, kernelHeight, kernelWidth);
			if (remainder) {
				MicrokernelRemainder(a, b, &result, row * kernelHeight, col * kernelWidth * 4, kernelHeight, kernelWidth);
			}
		}
	}
	// Do the remaining lines if any (we may still have some columns to fill in afterwards)
	for (size_t row = nRows * kernelHeight; row < a->rows; row++)
	{
		for (size_t col = 0; col < nCols; col++)
		{
			Microkernel(a, b, &result, row, col * kernelWidth * 4, 1, kernelWidth);
			if (remainder) {
				MicrokernelRemainder(a, b, &result, row * kernelHeight, col * kernelWidth * 4, 1, kernelWidth);
			}
		}
	}
	// Clean up what is left (columns)
	if (nCols * kernelWidth * 4 < b->cols)
	{
		for (size_t i = 0; i < a->rows; i++)
		{
			const float *aRow = a->contents + i * a->cols;
			float *resultRow = result.contents + i * result.cols;
			for (size_t k = 0; k < b->rows; k++)
			{
				const float *bRow = b->contents + k * b->cols;
				const float Aik = aRow[k];
				for (size_t j = nCols * kernelWidth * 4; j < b->cols; j++)
				{
					resultRow[j] += Aik * bRow[j];
				}
			}
		}
	}
	return result;
}

