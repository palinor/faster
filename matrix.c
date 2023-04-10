#pragma once
#ifdef __APPLE__
#include <arm_neon.h>
#else 
#include <immintrin.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "matrix.h"


void FreeMatrixf32(matrix_f32 *a) {
	free(a->contents);
}


matrix_f32 Matrixf32Copy(matrix_f32 *a) {
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = a->cols;
	result.ld = a->ld;
	result.contents = malloc(a->rows * a->cols * sizeof(float));
	if (!result.contents) {
		perror("Error allocating result.contents in Matrixf32Copy");
		exit(1);
	}
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
			float thisError = ABS(MatrixGetItem(a, i, j) - MatrixGetItem(b, i, j));
			if (!(thisError < MAX_FLOAT)) {
				return MAX_FLOAT; // handle NaN
			}
			maxError = (thisError > maxError) ? thisError : maxError;
		}
	}
	return maxError;
}

matrix_f32 Matrixf32Multiply(const matrix_f32 *a, const matrix_f32 *b)
{
	assert(a->cols == b->rows);
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = b->cols;
	result.contents = calloc(result.rows * result.cols, sizeof(float));
	if (!result.contents) {
		perror("Error allocating result.contents in Matrixf32Multiply");
		exit(1);
	}
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
	return result;
}

/*
Multiply inplace by a scalar
*/
void Matrixf32MultiplyInplace(matrix_f32 *a, float x) {
	for (size_t i = 0; i < a->rows; i++) {
		for (size_t j = 0; j < a->cols; j++) {
			MATRIX_ITEM(a, i, j) *= x;
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
	result.contents = malloc(a->rows * a->cols * sizeof(float));
	if (!result.contents) {
		perror("Error allocating results in Matrixf32Sum");
		exit(1);
	}
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
	if (!result.contents) {
		perror("Error allocating result.contents in RandomMatrixf32");
		exit(1);
	}
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			result.contents[i * cols + j] = ((float)rand()) / RAND_MAX;
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
	if (!result.contents) {
		perror("Error allocating result.contents in Onesf32");
		exit(1);
	}

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
void Microkernel(const matrix_f32 *a, const matrix_f32 *b, matrix_f32 *result, size_t resultRow, size_t resultCol, const size_t kernelWidth, const size_t kernelHeight)
{
	assert(a->cols == b->rows);
	__m256 resultElems[128];
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
				resultElems[i * kernelWidth + j] =  _mm256_fmadd_ps(broadcastA_i7, bElems_7j, resultElems[i * kernelWidth + j]);
			}
		}
	}

	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			_mm256_store_ps(MatrixGetAddr(result, resultRow + i, resultCol + SIMD_VECTOR_SIZE * j), resultElems[i * kernelWidth + j]);
		}
	}

}

/*
If k is not divisible by SIMD_VECTOR_SIZE, handle the addition of the remaining elements
*/

void MicrokernelRemainder(const matrix_f32 *a, const matrix_f32 *b, matrix_f32 *result, size_t resultRow, size_t resultCol, size_t kernelWidth, size_t kernelHeight)
{
	size_t remainder = a->cols % SIMD_VECTOR_SIZE;
	__m256 resultElems[SIMD_VECTOR_SIZE * 4];
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			resultElems[i * kernelWidth + j] = _mm256_load_ps(MatrixGetAddr(result, resultRow + i, resultCol + SIMD_VECTOR_SIZE * j));
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
			_mm256_store_ps(MatrixGetAddr(result, resultRow + i, resultCol + SIMD_VECTOR_SIZE * j), resultElems[i * kernelWidth + j]);
		}
	}
}
#endif
matrix_f32 Matrixf32MicrokernelMultiply(const matrix_f32 *a, const matrix_f32 *b, const size_t kernelWidth, const size_t kernelHeight)
{
	assert(a->cols == b->rows);
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = b->cols;
	result.contents = calloc(a->rows * b->cols, sizeof(float));
	const size_t nRows = a->rows / kernelHeight;
	const size_t nCols = b->cols / (kernelWidth * SIMD_VECTOR_SIZE);
	const size_t remainder = a->cols % SIMD_VECTOR_SIZE;
	// Do the main block
	for (size_t row = 0; row < nRows; row++)
	{
		for (size_t col = 0; col < nCols; col++)
		{
			Microkernel(a, b, &result, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, kernelHeight);
			if (remainder) {
				MicrokernelRemainder(a, b, &result, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, kernelHeight);
			}
		}
	}
	// Do the remaining lines if any (we may still have some columns to fill in afterwards)
	for (size_t row = nRows * kernelHeight; row < a->rows; row++)
	{
		for (size_t col = 0; col < nCols; col++)
		{
			Microkernel(a, b, &result, row, col * kernelWidth * SIMD_VECTOR_SIZE, 1, kernelWidth);
			if (remainder) {
				MicrokernelRemainder(a, b, &result, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, 1, kernelWidth);
			}
		}
	}
	// Clean up what is left (columns)
	if (nCols * kernelWidth * SIMD_VECTOR_SIZE < b->cols)
	{
		for (size_t i = 0; i < a->rows; i++)
		{
			const float *aRow = a->contents + i * a->cols;
			float *resultRow = result.contents + i * result.cols;
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
	return result;
}

