#pragma once

#ifdef _WIN32
#pragma warning( disable: 4305 4244 6386)
#endif
#include <stdio.h>

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

void Free(matrix_f32 *a);

matrix_f32 Matrixf32Copy(matrix_f32 *a);
float MatrixGetItem(matrix_f32 *a, size_t i, size_t j);
void MatrixPutItem(matrix_f32 *a, size_t i, size_t j, float x);
float *MatrixGetAddr(const matrix_f32 *a, size_t i, size_t j);
void PrintMatrixf32(matrix_f32 *result);
double GetMaxError(matrix_f32 *a, matrix_f32 *b);

matrix_f32 Matrixf32Multiply(const matrix_f32 *a, const matrix_f32 *b);
/*
Multiply inplace by a scalar
*/
void Matrixf32MultiplyInplace(matrix_f32 *a, float x);
matrix_f32 Matrixf32Sum(matrix_f32 *a, matrix_f32 *b);
matrix_f32 RandomMatrixf32(size_t rows, size_t cols);
matrix_f32 Onesf32(size_t rows, size_t cols);
matrix_f32 Identityf32(size_t n);
/*
Returns an upper triangular matrix with random non-zero elements
*/
matrix_f32 UpperRandomf32(size_t n);
/*
Returns an lower triangular matrix with random non-zero elements
*/
matrix_f32 LowerRandomf32(size_t n);
/*
If k is not divisible by SIMD_VECTOR_SIZE, handle the addition of the remaining elements
*/
matrix_f32 Matrixf32MicrokernelMultiply(const matrix_f32 *a, const matrix_f32 *b, const size_t kernelWidth, const size_t kernelHeight);
