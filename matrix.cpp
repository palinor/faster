#ifdef __APPLE__
#include <arm_neon.h>
#else 
#include <immintrin.h>
#include <omp.h>
#endif

#include <cassert>
#include <cstdlib>

#include <math.h>
#include "threadpool.cpp"

#define MAX_FLOAT 1e20
#define MIN_FLOAT 1e-12

struct Vectorf32 {
	size_t size = 0;
	float *contents = nullptr;
};

struct Vectorf64 {
	size_t size = 0;
	double *contents = nullptr;
};

struct Matrixf32
{
	size_t rows = 0;
	size_t cols = 0;
	float *contents = nullptr;
	size_t ld;
	float *operator[](int i) { return contents + i * ld; }
};

struct Matrixf64
{
	size_t rows = 0;
	size_t cols = 0;
	double *contents = nullptr;
	size_t ld;
	double *operator[](int i) { return contents + i * ld; }
};

struct MatrixTridiagonalf32 {
	size_t dimension = 0;
	float *upper_diagonal = nullptr;
	float *diagonal = nullptr;
	float *lower_diagonal = nullptr;
};

void MatrixTridiagonalf32MultiplyConstant(MatrixTridiagonalf32 *matrix, const float c) {
	matrix->upper_diagonal[0] *= c;
	matrix->diagonal[0] *= c;
	matrix->lower_diagonal[matrix->dimension - 1] *= c;
	matrix->diagonal[matrix->dimension - 1] *= c;
	for (size_t i = 1; i < matrix->dimension - 1; ++i) {
		matrix->upper_diagonal[i] *= c;
		matrix->diagonal[i] *= c;
		matrix->lower_diagonal[i] *= c;
	}
}

void MatrixTridiagonalAdd(MatrixTridiagonalf32 *output, MatrixTridiagonalf32 *a, MatrixTridiagonalf32 *b) {
	for (size_t i = 0; i < output->dimension; ++i) {
		output->upper_diagonal[i] = a->upper_diagonal[i] + b->upper_diagonal[i];
		output->diagonal[i] = a->diagonal[i] + b->diagonal[i];
		output->lower_diagonal[i] = a->lower_diagonal[i] + b->lower_diagonal[i];
	}
}


int MatrixTridiagonalf23ApplyToVector(Vectorf32 *output, MatrixTridiagonalf32 *matrix, Vectorf32 *vector) {
	if (matrix->dimension != vector->size) {
		return 1;
	}
	if (matrix->dimension != output->size) {
		return 2;
	}
	if (matrix->dimension < 3) {
		return 3;
	}
	if (output->contents == vector->contents) {
		return 4;
	}
	const size_t n = matrix->dimension;
	output->contents[0] = matrix->upper_diagonal[0] * vector->contents[1] + matrix->diagonal[0] * vector->contents[0];
	output->contents[n - 1] = matrix->lower_diagonal[n - 1] * vector->contents[n - 2] + matrix->diagonal[n - 1] * vector->contents[n - 1];
	for (size_t i = 1; i < n - 1; ++i) {
		output->contents[i] = matrix->upper_diagonal[i] * vector->contents[i + 1]
			+ matrix->diagonal[i] * vector->contents[i]
			+ matrix->lower_diagonal[i] * vector->contents[i - 1];
	}
	return 0;
}

int MatrixTridiagonalf32InvertEquation(Vectorf32 *output, MatrixTridiagonalf32 *matrix, Vectorf32 *vector, Vectorf32 *scratch) {
	if (matrix->dimension != vector->size) {
		return 1;
	}
	if (matrix->dimension != output->size) {
		return 2;
	}
	if (matrix->dimension < 3) {
		return 3;
	}
	if (output->contents == vector->contents) {
		return 4;
	}

	const size_t n = matrix->dimension;

	scratch->contents[0] = matrix->upper_diagonal[0] / matrix->diagonal[0];
	vector->contents[0] = vector->contents[0] / matrix->diagonal[0];
	for (size_t i = 1; i < matrix->dimension; ++i) {
		scratch->contents[i] = matrix->upper_diagonal[i] / (matrix->diagonal[i] - matrix->lower_diagonal[i] * scratch->contents[i - 1]);
		vector->contents[i] -= matrix->lower_diagonal[i] * vector->contents[i - 1];
		vector->contents[i] /= matrix->diagonal[i] - matrix->lower_diagonal[i] * scratch->contents[i - 1];
	}

	output->contents[n - 1] = vector->contents[n - 1];
	for (int i = n - 2; i > 0; --i) {
		output->contents[i] = vector->contents[i] - scratch->contents[i] * output->contents[i + 1];
	}
	return 0;
}


inline void matrixf32Free(Matrixf32 *a) {
	free(a->contents);
}

inline void matrixf64Free(Matrixf64 *a) {
	free(a->contents);
}

inline float matrixf32GetItem(Matrixf32 *a, size_t i, size_t j)
{
	return a->contents[i * a->ld + j];
}

inline double matrixf64GetItem(Matrixf64 *a, size_t i, size_t j)
{
	return a->contents[i * a->ld + j];
}

inline void matrixf32PutItem(Matrixf32 *a, size_t i, size_t j, float x) {
	a->contents[i * a->ld + j] = x;
}

inline void matrixf64PutItem(Matrixf64 *a, size_t i, size_t j, float x) {
	a->contents[i * a->ld + j] = x;
}

inline float *matrixf32GetAddr(const Matrixf32 *a, size_t i, size_t j)
{
	return a->contents + (i * a->ld + j);
}

inline double *matrixf64GetAddr(const Matrixf64 *a, size_t i, size_t j)
{
	return a->contents + (i * a->ld + j);
}

inline float f32Abs(float x) {
	if (x < 0) {
		return -x;
	}
	else {
		return x;
	}
}

inline double f64Abs(double x) {
	if (x < 0) {
		return -x;
	}
	else {
		return x;
	}
}

inline int intAbs(int x) {
	if (x < 0) {
		return -x;
	}
	else {
		return x;
	}
}

// todo(AION) test this out (step through with debugger to make sure it does what we expect)
int matrixf32Copy(Matrixf32 *result, Matrixf32 *input) {
	if (!result->contents) {
		result->contents = (float *)malloc(input->rows * input->cols * sizeof(float));
		if (!result->contents) return 1;
		result->rows = input->rows;
		result->cols = input->cols;
		result->ld = input->ld;
	}
	float *result_element = result->contents;
	float *input_element = input->contents;
	for (size_t i = 0; i < input->rows; ++i) {
		for (size_t j = 0; j < input->cols; ++j) {
			*result_element++ = *input_element++;
		}
	}
	return 0;
}

int matrixf64Copy(Matrixf64 *result, Matrixf64 *input) {
	if (!result->contents) {
		result->contents = (double *)malloc(input->rows * input->cols * sizeof(double));
		if (!result->contents) return 1;
		result->rows = input->rows;
		result->cols = input->cols;
		result->ld = input->ld;
	}
	double *result_element = result->contents;
	double *input_element = input->contents;
	for (size_t i = 0; i < input->rows; ++i) {
		for (size_t j = 0; j < input->cols; ++j) {
			*result_element++ = *input_element++;
		}
	}
	return 0;
}


void matrixf32Print(Matrixf32 *result)
{
	for (size_t i = 0; i < result->rows; ++i)
	{
		for (size_t j = 0; j < result->cols; ++j)
		{
			printf("%f ", matrixf32GetItem(result, i, j));
		}
		printf("\n");
	}
}

void matrixf64Print(Matrixf64 *result)
{
	for (size_t i = 0; i < result->rows; ++i)
	{
		for (size_t j = 0; j < result->cols; ++j)
		{
			printf("%f ", matrixf64GetItem(result, i, j));
		}
		printf("\n");
	}
}

float matrixf32Max(Matrixf32 *a)
{
	float max = 0;
	float *a_element = a->contents;
	for (size_t i = 0; i < a->rows; ++i)
	{
		for (size_t j = 0; j < a->cols; ++j)
		{
			float abs_value = f32Abs(*a_element++);
			if (!(abs_value < MAX_FLOAT)) {
				return (float)MAX_FLOAT; // handle NaN
			}
			max = (abs_value > max) ? abs_value : max;
		}
	}
	return max;
}

double matrixf64Max(Matrixf64 *a)
{
	double max = 0;
	double *a_element = a->contents;
	for (size_t i = 0; i < a->rows; ++i)
	{
		for (size_t j = 0; j < a->cols; ++j)
		{
			double abs_value = f64Abs(*a_element++);
			if (!(abs_value < MAX_FLOAT)) {
				return MAX_FLOAT; // handle NaN
			}
			max = (abs_value > max) ? abs_value : max;
		}
	}
	return max;
}


int matrixf32MultiplyToTarget(Matrixf32 *result, const Matrixf32 *a, const Matrixf32 *b)
{
	// check multiplication is allowed
	if (a->cols != b->rows) return 1;
	// if result has not been allocated yet, we try and allocate ourselves
	if (!result->contents) {
		result->contents = (float *)malloc(a->rows * b->cols * sizeof(float));
		if (!result->contents) return 1;
		result->rows = a->rows;
		result->cols = b->cols;
		result->ld = result->cols;
	}
	float *a_row = a->contents;
	float *result_row = result->contents;
	for (size_t i = 0; i < a->rows; ++i)
	{
		float *b_row = b->contents;
		float *A_ik = a_row;
		for (size_t k = 0; k < b->rows; ++k)
		{
			// #pragma clang loop vectorize_width(8) // interleave_count(4)
			for (size_t j = 0; j < result->cols; j++)
			{
#ifdef __APPLE__
				result_row[j] += *A_ik * b_row[j];
				/*
				This part is actually out of line with the microkernel loop (no difference in -Ofast, but visible in other compile modes)
				The above lines up just fine in all

				BUUUUUUUT.... not on windows. So we actually need to swap depending on the compiler. Lol.
				Maybe it was the conversion float->double->float? If Apple math header has a float version of fma as opposed to fmaf
				*/
#else
				result_row[j] = fmaf(b_row[j], *A_ik, result_row[j]);
#endif
			}
			b_row += b->ld;
			++A_ik;
		}
		a_row += a->ld;
		result_row += result->ld;
	}
	return 0;
}

int matrixf64MultiplyToTarget(Matrixf64 *result, const Matrixf64 *a, const Matrixf64 *b)
{
	// Check multiplication is allowed
	if (a->cols != b->rows) return 1;
	// if the result has not been allocated, we try and allocate here
	if (!result->contents) {
		result->contents = (double *)malloc(a->rows * b->cols * sizeof(double));
		if (!result->contents) return 1;
		result->rows = a->rows;
		result->cols = b->cols;
		result->ld = result->cols;
	}
	double *a_row = a->contents;
	double *result_row = result->contents;
	for (size_t i = 0; i < a->rows; ++i)
	{
		double *b_row = b->contents;
		double *A_ik = a_row;
		for (size_t k = 0; k < b->rows; k++)
		{
			// #pragma clang loop vectorize_width(8) // interleave_count(4)
			for (size_t j = 0; j < result->cols; ++j)
			{
#ifdef __APPLE__
				result_row[j] += *A_ik * b_row[j];
				/*
				This part is actually out of line with the microkernel loop (no difference in -Ofast, but visible in other compile modes)
				The above lines up just fine in all

				BUUUUUUUT.... not on windows. So we actually need to swap depending on the compiler. Lol.
				*/
#else
				result_row[j] = fma(b_row[j], *A_ik, result_row[j]);
#endif
			}
			b_row += b->ld;
			++A_ik;
		}
		a_row += a->ld;
		result_row += result->ld;
	}
	return 0;
}

Matrixf32 matrixf32Multiply(const Matrixf32 *a, const Matrixf32 *b)
{
	assert(a->cols == b->rows);
	Matrixf32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = b->cols;
	result.contents = (float *)calloc(result.rows * result.cols, sizeof(float));
	matrixf32MultiplyToTarget(&result, a, b);
	return result;
}

Matrixf64 matrixf64Multiply(const Matrixf64 *a, const Matrixf64 *b)
{
	assert(a->cols == b->rows);
	Matrixf64 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = b->cols;
	result.contents = (double *)calloc(result.rows * result.cols, sizeof(double));
	matrixf64MultiplyToTarget(&result, a, b);
	return result;
}

/*
Multiply inplace by a scalar
*/
void matrixf32MultiplyInplace(Matrixf32 *a, float x) {
	float *a_element = a->contents;
	for (size_t i = 0; i < a->rows; ++i) {
		for (size_t j = 0; j < a->cols; ++j) {
			*a_element++ *= x;
		}
	}
}

void matrixf64MultiplyInplace(Matrixf64 *a, double x) {
	double *a_element = a->contents;
	for (size_t i = 0; i < a->rows; i++) {
		for (size_t j = 0; j < a->cols; j++) {
			*a_element++ *= x;
		}
	}
}

int matrixf32Sum(Matrixf32 *result, Matrixf32 *a, Matrixf32 *b) {
	if (a->rows != b->rows || a->cols != b->cols) return 1;
	if (!result->contents) {
		result->contents = (float *)malloc(a->rows * a->cols * sizeof(float));
		if (!result->contents) return 1;
		result->rows = a->rows;
		result->cols = a->cols;
		result->ld = result->cols;
	}
	float *result_element = result->contents;
	float *a_element = a->contents;
	float *b_element = b->contents;
	for (size_t i = 0; i < result->rows; i++) {
		for (size_t j = 0; j < result->cols; j++) {
			*result_element++ = *a_element++ + *b_element++;
		}
	}
	return 0;
}

int matrixf64Sum(Matrixf64 *result, Matrixf64 *a, Matrixf64 *b) {
	if (a->rows != b->rows || a->cols != b->cols) return 1;
	if (!result->contents) {
		result->contents = (double *)malloc(a->rows * a->cols * sizeof(double));
		if (!result->contents) return 1;
		result->rows = a->rows;
		result->cols = a->cols;
		result->ld = result->cols;
	}
	double *result_element = result->contents;
	double *a_element = a->contents;
	double *b_element = b->contents;
	for (size_t i = 0; i < result->rows; i++) {
		for (size_t j = 0; j < result->cols; j++) {
			*result_element++ = *a_element++ + *b_element++;
		}
	}
	return 0;
}

int matrixf32RandomInit(Matrixf32 *result, size_t rows, size_t cols) {
	// let's make sure this memory does not exist before we re-allocate it
	if (result->contents) free(result->contents);
	result->contents = (float *)malloc(rows * cols * sizeof(float));
	if (!result->contents) return 1;
	result->rows = rows;
	result->cols = cols;
	result->ld = cols;
	float *result_element = result->contents;
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < cols; ++j)
		{
			int random_int = rand();
			float random_float = ((float)random_int) / RAND_MAX;
			*result_element++ = random_float;
		}
	}
	return 0;
}


Matrixf32 matrixf32Random(size_t rows, size_t cols)
{
	Matrixf32 result;
	int err = matrixf32RandomInit(&result, rows, cols);
	if (err) {
		perror("Error allocating matrixf32");
		exit(1);
	}
	return result;
}

int matrixf64RandomInit(Matrixf64 *result, size_t rows, size_t cols) {
	// Let's make sure this does not exist
	if (result->contents) free(result->contents);
	result->contents = (double *)malloc(rows * cols * sizeof(double));
	if (!result->contents) return 1;
	result->rows = rows;
	result->cols = cols;
	result->ld = cols;
	double *result_element = result->contents;
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < cols; ++j)
		{
			int random_int = rand();
			double random_float = ((double)random_int) / RAND_MAX;
			*result_element++ = random_float;
		}
	}
	return 0;
}

Matrixf64 matrixf64Random(size_t rows, size_t cols)
{
	Matrixf64 result;
	int err = matrixf64RandomInit(&result, rows, cols);
	if (err) {
		perror("Error allocating matrixf64");
		exit(1);
	}
	return result;
}

int matrixf32InitOnes(Matrixf32 *result, size_t rows, size_t cols) {
	if (result->contents) free(result->contents);
	result->contents = (float *)malloc(rows * cols * sizeof(float));
	if (!result->contents) return 1;
	result->rows = rows;
	result->cols = cols;
	result->ld = cols;
	float *result_element = result->contents;
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			*result_element++ = 1;
		}
	}
	return 0;
}

Matrixf32 matrixf32Ones(size_t rows, size_t cols)
{
	Matrixf32 result;
	int err = matrixf32InitOnes(&result, rows, cols);
	if (err) {
		perror("Error allocating matrixf32");
		exit(1);
	}
	return result;
}

int matrixf64InitOnes(Matrixf64 *result, size_t rows, size_t cols) {
	if (result->contents) free(result->contents);
	result->contents = (double *)malloc(rows * cols * sizeof(double));
	if (!result->contents) return 1;
	result->rows = rows;
	result->cols = cols;
	result->ld = cols;
	double *result_element = result->contents;
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			*result_element++ = 1;
		}
	}
	return 0;
}

Matrixf64 matrixf64Ones(size_t rows, size_t cols)
{
	Matrixf64 result;
	int err = matrixf64InitOnes(&result, rows, cols);
	if (err) {
		perror("Error allocating matrixf32");
		exit(1);
	}
	return result;
}

int matrixf32InitIdentity(Matrixf32 *result, size_t n) {
	if (result->contents) free(result->contents);
	result->contents = (float *)calloc(n * n, sizeof(float));
	if (!result->contents) return 1;
	result->rows = n;
	result->cols = n;
	result->ld = n;
	float *result_element = result->contents;
	for (size_t i = 0; i < n; ++i) {
		*result_element = 1;
		result_element += n + 1;
	}
	return 0;
}

Matrixf32 Matrixf32CreateIdentity(size_t n)
{
	Matrixf32 result;
	result.rows = n;
	result.cols = n;
	result.ld = n;
	result.contents = (float *)calloc(n * n, sizeof(float));
	float *result_element = result.contents;
	for (size_t i = 0; i < n; i++)
	{
		*result_element = 1;
		result_element += n + 1;
	}
	return result;
}

int matrixf32InitUpperRandom(Matrixf32 *result, size_t n) {
	if (result->contents) free(result->contents);
	result->contents = (float *)calloc(n * n, sizeof(float));
	if (!result->contents) return 1;
	result->rows = n;
	result->cols = n;
	result->ld = n;
	float *result_element = result->contents;
	for (size_t i = 0; i < n; ++i) {
		result_element += i;
		for (size_t j = i; j < n; ++j) {
			int random_int = rand();
			float random_float = ((float)random_int) / RAND_MAX;
			*result_element++ = random_float;
		}
	}
	return 0;
}

/*
Returns an upper triangular matrix with random non-zero elements
*/
Matrixf32 matrixf32CreateUpperRandom(size_t n) {
	Matrixf32 result;
	result.rows = n;
	result.cols = n;
	result.ld = n;
	int err = matrixf32InitUpperRandom(&result, n);
	if (err) {
		perror("Failure in allocating Matrixf32");
		exit(EXIT_FAILURE);
	}
	return result;
}

int matrixf32InitLowerRandom(Matrixf32 *result, size_t n) {
	if (result->contents) free(result->contents);
	result->contents = (float *)calloc(n * n, sizeof(float));
	if (!result->contents) return 1;
	result->rows = n;
	result->cols = n;
	result->ld = n;
	float *result_element = result->contents;
	for (size_t i = 0; i < n; ++i) {
		result_element += i;
		for (size_t j = i + 1; j > 0; --j) {
			int random_int = rand();
			float random_float = ((float)random_int) / RAND_MAX;
			*result_element++ = random_float;
		}
	}
	return 0;
}

int matrixf64InitLowerRandom(Matrixf64 *result, size_t n) {
	if (result->contents) free(result->contents);
	result->contents = (double *)calloc(n * n, sizeof(double));
	if (!result->contents) return 1;
	result->rows = n;
	result->cols = n;
	result->ld = n;
	double *result_element = result->contents;
	for (size_t i = 0; i < n; ++i) {
		result_element += i;
		for (size_t j = i + 1; j > 0; --j) {
			int random_int = rand();
			double random_float = ((double)random_int) / RAND_MAX;
			*result_element++ = random_float;
		}
	}
	return 0;
}
/*
Returns an lower triangular matrix with random non-zero elements
*/
Matrixf32 LowerRandomf32(size_t n) {
	Matrixf32 result;
	int err = matrixf32InitLowerRandom(&result, n);
	if (err) {
		perror("Allocating Matrixf32 failed");
		exit(1);
	}
	return result;
}

/*
Internal block for matrix multiplication - compute the result in blocks of isze kernel_height x (4 * kernel_width)
*/
#ifdef __APPLE__
// ARM implimentation
#define SIMD_VECTOR_SIZE 4
int microkernel(Matrixf32 *result_matrix, const Matrixf32 *a, const Matrixf32 *b, size_t result_row, size_t result_col, size_t kernel_width, size_t kernel_height)
{
	if (a->cols != b->rows) return 1;
	if (kernel_width * kernel_height < 128) return 1;
	float32x4_t result_elements[128];
	for (size_t i = 0; i < kernel_height; i++)
	{
		for (size_t j = 0; j < kernel_width; j++)
		{
			result_elements[i * kernel_width + j] = vdupq_n_f32(0);
		}
	}
	for (size_t k = 0; k < a->cols / SIMD_VECTOR_SIZE; k++)
	{
		for (size_t i = 0; i < kernel_height; i++)
		{
			float32x4_t a_element_i = vld1q_f32(matrixf32GetAddr(a, result_row + i, SIMD_VECTOR_SIZE * k));
			for (size_t j = 0; j < kernel_width; j++)
			{
				float32x4_t b_element_i0 = vld1q_f32(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k, result_col + SIMD_VECTOR_SIZE * j));
				result_elements[i * kernel_width + j] = vmlaq_laneq_f32(result_elements[i * kernel_width + j], b_element_i0, a_element_i, 0);

				float32x4_t b_element_i1 = vld1q_f32(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k + 1, result_col + SIMD_VECTOR_SIZE * j));
				result_elements[i * kernel_width + j] = vmlaq_laneq_f32(result_elements[i * kernel_width + j], b_element_i1, a_element_i, 1);

				float32x4_t b_element_i2 = vld1q_f32(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k + 2, result_col + SIMD_VECTOR_SIZE * j));
				result_elements[i * kernel_width + j] = vmlaq_laneq_f32(result_elements[i * kernel_width + j], b_element_i2, a_element_i, 2);

				float32x4_t b_element_i3 = vld1q_f32(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k + 3, result_col + SIMD_VECTOR_SIZE * j));
				result_elements[i * kernel_width + j] = vmlaq_laneq_f32(result_elements[i * kernel_width + j], b_element_i3, a_element_i, 3);
			}
		}
	}
	// If k is not divisible by SIMD_VECTOR_SIZE, handle the remainder here
	size_t remainder = a->cols % SIMD_VECTOR_SIZE;
	if (remainder) {
		for (size_t i = 0; i < kernel_height; i++) {
			float32x4_t a_remainder_element_i = vld1q_f32(matrixf32GetAddr(a, result_row + i, a->cols - remainder));
			for (size_t j = 0; j < kernel_width; j++) {
				float32x4_t b_remainder_element_0j = vld1q_f32(matrixf32GetAddr(b, a->cols - remainder, result_col + SIMD_VECTOR_SIZE * j));
				result_elements[i * kernel_width + j] = vmlaq_laneq_f32(result_elements[i * kernel_width + j], b_remainder_element_0j, a_remainder_element_i, 0);
				if (remainder > 1) {
					float32x4_t b_remainder_element_1j = vld1q_f32(matrixf32GetAddr(b, a->cols - remainder + 1, result_col + SIMD_VECTOR_SIZE * j));
					result_elements[i * kernel_width + j] = vmlaq_laneq_f32(result_elements[i * kernel_width + j], b_remainder_element_1j, a_remainder_element_i, 1);
				}
				if (remainder > 2) {
					float32x4_t b_remainder_element_2j = vld1q_f32(matrixf32GetAddr(b, a->cols - remainder + 2, result_col + SIMD_VECTOR_SIZE * j));
					result_elements[i * kernel_width + j] = vmlaq_laneq_f32(result_elements[i * kernel_width + j], b_remainder_element_2j, a_remainder_element_i, 2);
				}
			}
		}
	}
	float *result_location = matrixf32GetAddr(result_matrix, result_row, result_col);
	for (size_t i = 0; i < kernel_height; i++)
	{
		for (size_t j = 0; j < kernel_width; j++)
		{
			vst1q_f32(result_location + i * b->cols + j * SIMD_VECTOR_SIZE, result_elements[i * kernel_width + j]);
		}
	}
	return 0;
}

/*
If k is not divisible by SIMD_VECTOR_SIZE, handle the addition of the remaining elements

void MicrokernelRemainder(const Matrixf32 *a, const Matrixf32 *b, float *result_location, size_t result_row, size_t result_column, size_t kernel_width, size_t kernel_height)
{
	size_t remainder = a->cols % SIMD_VECTOR_SIZE;
	float32x4_t result_elements[32];
	for (size_t i = 0; i < kernel_height; i++)
	{
		for (size_t j = 0; j < kernel_width; j++)
		{
			result_elements[i * kernel_width + j] = vdupq_n_f32(0);
		}
	}
	for (size_t i = 0; i < kernel_height; i++)
	{
		float32x4_t a_remainder_element_i = vld1q_f32(matrixf32GetAddr(a, result_row + i, a->cols - remainder));
		for (size_t j = 0; j < kernel_width; j++)
		{
			float32x4_t b_remainder_element_0j = vld1q_f32(matrixf32GetAddr(b, a->cols - remainder, result_column + SIMD_VECTOR_SIZE * j));
			result_elements[i * kernel_width + j] = vmlaq_laneq_f32(result_elements[i * kernel_width + j], b_remainder_element_0j, a_remainder_element_i, 0);
			if (remainder > 1)
			{
				float32x4_t b_remainder_element_1j = vld1q_f32(matrixf32GetAddr(b, a->cols - remainder + 1, result_column + SIMD_VECTOR_SIZE * j));
				result_elements[i * kernel_width + j] = vmlaq_laneq_f32(result_elements[i * kernel_width + j], b_remainder_element_1j, a_remainder_element_i, 1);
			}
			if (remainder > 2)
			{
				float32x4_t b_remainder_element_2j = vld1q_f32(matrixf32GetAddr(b, a->cols - remainder + 2, result_column + SIMD_VECTOR_SIZE * j));
				result_elements[i * kernel_width + j] = vmlaq_laneq_f32(result_elements[i * kernel_width + j], b_remainder_element_2j, a_remainder_element_i, 2);
			}
		}
	}
	for (size_t i = 0; i < kernel_height; i++)
	{
		for (size_t j = 0; j < kernel_width; j++)
		{
			vst1q_f32(result_location + i * b->cols + j * SIMD_VECTOR_SIZE, result_elements[i * kernel_width + j]);
		}
	}
} */
#else
// AMD64 implementation
#define SIMD_VECTOR_SIZE 8 
void matrixf32MicrokernelApply(Matrixf32 *result_matrix, const Matrixf32 *a, const Matrixf32 *b, size_t result_row, size_t result_column, const size_t kernel_width, const size_t kernel_height)
{
	assert(a->cols == b->rows);
	assert(kernel_height * kernel_width < 32);
	__m256 result_elements[32];
	for (size_t i = 0; i < kernel_height; i++)
	{
		for (size_t j = 0; j < kernel_width; j++)
		{
			result_elements[i * kernel_width + j] = _mm256_set1_ps(0);
		}
	}
	for (size_t k = 0; k < a->cols / SIMD_VECTOR_SIZE; k++)
	{
		for (size_t i = 0; i < kernel_height; i++)
		{
			float *aElems_i = matrixf32GetAddr(a, result_row + i, SIMD_VECTOR_SIZE * k);
			for (size_t j = 0; j < kernel_width; j++)
			{
				//todo(AION): we have 3 fma ports, why can't we use 3 accumulators? It seemed needlessly slow when I tried 2
				__m256 b_elements_0j = _mm256_loadu_ps(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k, result_column + SIMD_VECTOR_SIZE * j));
				__m256 a_element_broadcast_i0 = _mm256_set1_ps(*aElems_i);
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(a_element_broadcast_i0, b_elements_0j, result_elements[i * kernel_width + j]);

				__m256 b_elements_1j = _mm256_loadu_ps(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k + 1, result_column + SIMD_VECTOR_SIZE * j));
				__m256 a_element_broadcast_i1 = _mm256_set1_ps(*(aElems_i + 1));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(a_element_broadcast_i1, b_elements_1j, result_elements[i * kernel_width + j]);

				__m256 b_elements_2j = _mm256_loadu_ps(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k + 2, result_column + SIMD_VECTOR_SIZE * j));
				__m256 a_element_broadcast_i2 = _mm256_set1_ps(*(aElems_i + 2));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(a_element_broadcast_i2, b_elements_2j, result_elements[i * kernel_width + j]);

				__m256 b_elements_3j = _mm256_loadu_ps(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k + 3, result_column + SIMD_VECTOR_SIZE * j));
				__m256 a_element_broadcast_i3 = _mm256_set1_ps(*(aElems_i + 3));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(a_element_broadcast_i3, b_elements_3j, result_elements[i * kernel_width + j]);

				__m256 b_elements_4j = _mm256_loadu_ps(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k + 4, result_column + SIMD_VECTOR_SIZE * j));
				__m256 a_element_broadcast_i4 = _mm256_set1_ps(*(aElems_i + 4));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(a_element_broadcast_i4, b_elements_4j, result_elements[i * kernel_width + j]);

				__m256 b_elements_5j = _mm256_loadu_ps(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k + 5, result_column + SIMD_VECTOR_SIZE * j));
				__m256 a_element_broadcast_i5 = _mm256_set1_ps(*(aElems_i + 5));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(a_element_broadcast_i5, b_elements_5j, result_elements[i * kernel_width + j]);

				__m256 b_elements_6j = _mm256_loadu_ps(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k + 6, result_column + SIMD_VECTOR_SIZE * j));
				__m256 a_element_broadcast_i6 = _mm256_set1_ps(*(aElems_i + 6));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(a_element_broadcast_i6, b_elements_6j, result_elements[i * kernel_width + j]);

				__m256 b_elements_7j = _mm256_loadu_ps(matrixf32GetAddr(b, SIMD_VECTOR_SIZE * k + 7, result_column + SIMD_VECTOR_SIZE * j));
				__m256 a_element_broadcast_i7 = _mm256_set1_ps(*(aElems_i + 7));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(a_element_broadcast_i7, b_elements_7j, result_elements[i * kernel_width + j]);
			}
		}
	}
	// If k is not divisible by SIMD_VECTOR_SIZE, handle the remainder
	size_t remainder = a->cols % SIMD_VECTOR_SIZE;
	if (remainder) {
		for (size_t i = 0; i < kernel_height; i++)
		{
			float *a_remainder_element_i = matrixf32GetAddr(a, result_row + i, a->cols - remainder);
			for (size_t j = 0; j < kernel_width; j++)
			{
				__m256 b_remainder_element_0j = _mm256_loadu_ps(matrixf32GetAddr(b, a->cols - remainder, result_column + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i0 = _mm256_set1_ps(*a_remainder_element_i);
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i0, b_remainder_element_0j, result_elements[i * kernel_width + j]);
				if (remainder > 1)
				{
					__m256 b_remainder_element_1j = _mm256_loadu_ps(matrixf32GetAddr(b, a->cols - remainder + 1, result_column + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i1 = _mm256_set1_ps(*(a_remainder_element_i + 1));
					result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i1, b_remainder_element_1j, result_elements[i * kernel_width + j]);
				}
				if (remainder > 2)
				{
					__m256 b_remainder_element_2j = _mm256_loadu_ps(matrixf32GetAddr(b, a->cols - remainder + 2, result_column + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i2 = _mm256_set1_ps(*(a_remainder_element_i + 2));
					result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i2, b_remainder_element_2j, result_elements[i * kernel_width + j]);
				}
				if (remainder > 3)
				{
					__m256 b_remainder_element_3j = _mm256_loadu_ps(matrixf32GetAddr(b, a->cols - remainder + 3, result_column + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i3 = _mm256_set1_ps(*(a_remainder_element_i + 3));
					result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i3, b_remainder_element_3j, result_elements[i * kernel_width + j]);
				}
				if (remainder > 4)
				{
					__m256 b_remainder_element_4j = _mm256_loadu_ps(matrixf32GetAddr(b, a->cols - remainder + 4, result_column + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i4 = _mm256_set1_ps(*(a_remainder_element_i + 4));
					result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i4, b_remainder_element_4j, result_elements[i * kernel_width + j]);
				}
				if (remainder > 5)
				{
					__m256 b_remainder_element_5j = _mm256_loadu_ps(matrixf32GetAddr(b, a->cols - remainder + 5, result_column + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i5 = _mm256_set1_ps(*(a_remainder_element_i + 5));
					result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i5, b_remainder_element_5j, result_elements[i * kernel_width + j]);
				}
				if (remainder > 6)
				{
					__m256 b_remainder_element_6j = _mm256_loadu_ps(matrixf32GetAddr(b, a->cols - remainder + 6, result_column + SIMD_VECTOR_SIZE * j));
					__m256 aRemainderElems_i6 = _mm256_set1_ps(*(a_remainder_element_i + 6));
					result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i6, b_remainder_element_6j, result_elements[i * kernel_width + j]);
				}
			}
		}
	}
	float *result_location = matrixf32GetAddr(result_matrix, result_row, result_column);
	for (size_t i = 0; i < kernel_height; i++)
	{
		for (size_t j = 0; j < kernel_width; j++)
		{
			_mm256_storeu_ps(result_location + i * b->cols + j * SIMD_VECTOR_SIZE, result_elements[i * kernel_width + j]);
		}
	}

}

/*
If k is not divisible by SIMD_VECTOR_SIZE, handle the addition of the remaining elements
void MicrokernelRemainder(const Matrixf32 *a, const Matrixf32 *b, float *result_location, size_t result_row, size_t result_column, size_t kernel_width, size_t kernel_height)
{
	size_t remainder = a->cols % SIMD_VECTOR_SIZE;
	__m256 result_elements[32]:
	for (size_t i = 0; i < kernel_height; i++)
	{
		for (size_t j = 0; j < kernel_width; j++)
		{
			result_elements[i * kernel_width + j] = _mm256_loadu_ps(0);
		}
	}
	for (size_t i = 0; i < kernel_height; i++)
	{
		float *a_remainder_element_i = matrixGetAddr(a, result_row + i, a->cols - remainder);
		for (size_t j = 0; j < kernel_width; j++)
		{
			__m256 b_remainder_element_0j = _mm256_loadu_ps(matrixGetAddr(b, a->cols - remainder, result_column + SIMD_VECTOR_SIZE * j));
			__m256 aRemainderElems_i0 = _mm256_set1_ps(*a_remainder_element_i);
			result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i0, b_remainder_element_0j, result_elements[i * kernel_width + j]);
			if (remainder > 1)
			{
				__m256 b_remainder_element_1j = _mm256_loadu_ps(matrixGetAddr(b, a->cols - remainder + 1, result_column + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i1 = _mm256_set1_ps(*(a_remainder_element_i + 1));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i1, b_remainder_element_1j, result_elements[i * kernel_width + j]);
			}
			if (remainder > 2)
			{
				__m256 b_remainder_element_2j = _mm256_loadu_ps(matrixGetAddr(b, a->cols - remainder + 2, result_column + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i2 = _mm256_set1_ps(*(a_remainder_element_i + 2));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i2, b_remainder_element_2j, result_elements[i * kernel_width + j]);
			}
			if (remainder > 3)
			{
				__m256 b_remainder_element_3j = _mm256_loadu_ps(matrixGetAddr(b, a->cols - remainder + 3, result_column + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i3 = _mm256_set1_ps(*(a_remainder_element_i + 3));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i3, b_remainder_element_3j, result_elements[i * kernel_width + j]);
			}
			if (remainder > 4)
			{
				__m256 b_remainder_element_4j = _mm256_loadu_ps(matrixGetAddr(b, a->cols - remainder + 4, result_column + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i4 = _mm256_set1_ps(*(a_remainder_element_i + 4));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i4, b_remainder_element_4j, result_elements[i * kernel_width + j]);
			}
			if (remainder > 5)
			{
				__m256 b_remainder_element_5j = _mm256_loadu_ps(matrixGetAddr(b, a->cols - remainder + 5, result_column + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i5 = _mm256_set1_ps(*(a_remainder_element_i + 5));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i5, b_remainder_element_5j, result_elements[i * kernel_width + j]);
			}
			if (remainder > 6)
			{
				__m256 b_remainder_element_6j = _mm256_loadu_ps(matrixGetAddr(b, a->cols - remainder + 6, result_column + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i6 = _mm256_set1_ps(*(a_remainder_element_i + 6));
				result_elements[i * kernel_width + j] = _mm256_fmadd_ps(aRemainderElems_i6, b_remainder_element_6j, result_elements[i * kernel_width + j]);
			}
		}
	}
	for (size_t i = 0; i < kernel_height; i++)
	{
		for (size_t j = 0; j < kernel_width; j++)
		{
			_mm256_store_ps(result_location, result_elements[i * kernel_width + j]);
		}
	}
}*/
#endif

struct Matrixf32MicrokernelTaskInfo {
	const ConcurrentTaskQueue *queue;
	const Matrixf32 *left_matrix;
	const Matrixf32 *right_matrix;
	Matrixf32 *result_matrix;
	size_t result_row_number;
	size_t kernel_width;
	size_t kernel_height;
	size_t task_number;
};


void matrixMicrokernelRowTask(void *microkernel_task_infoInput) {
	Matrixf32MicrokernelTaskInfo *microkernel_task_info = reinterpret_cast<Matrixf32MicrokernelTaskInfo *>(microkernel_task_infoInput);
	const Matrixf32 *a = microkernel_task_info->left_matrix;
	const Matrixf32 *b = microkernel_task_info->right_matrix;
	size_t row = microkernel_task_info->result_row_number;
	size_t kernel_width = microkernel_task_info->kernel_width;
	size_t kernel_height = microkernel_task_info->kernel_height;
	Matrixf32 *result_matrix = microkernel_task_info->result_matrix;
	size_t number_of_columns = microkernel_task_info->right_matrix->cols / (kernel_width * SIMD_VECTOR_SIZE);
	for (size_t col = 0; col < number_of_columns; col++) {
		matrixf32MicrokernelApply(result_matrix, a, b, row * kernel_height, col * kernel_width * SIMD_VECTOR_SIZE, kernel_width, kernel_height);
	}
	if (microkernel_task_info->queue){
		microkernel_task_info->queue->task_is_done[microkernel_task_info->task_number] = true;
	}
}



void Matrixf32MicrokernelMultiply(
	Matrixf32 *result,
	const Matrixf32 *a,
	const Matrixf32 *b,
	const size_t kernel_width,
	const size_t kernel_height,
	ThreadPool *pool
)
{
	assert(a->cols == b->rows);
	assert(result->rows == a->rows);
	assert(result->cols == b->cols);
	assert(!!result->contents);
	const size_t number_of_rows = a->rows / kernel_height;
	const size_t remaining_rows_after_microkernel = a->rows % kernel_height;
	const size_t number_of_columns = b->cols / (kernel_width * SIMD_VECTOR_SIZE);
	bool skip_main_block = number_of_columns == 0; // if we don't have enough columns, skip the vectorized block
	//TaskHandle *futures = (TaskHandle *)malloc(sizeof(TaskHandle) * number_of_rows);
	// Do the main block
	if (!skip_main_block) {
		Matrixf32MicrokernelTaskInfo *microkernel_task_info_list = (Matrixf32MicrokernelTaskInfo *)malloc(sizeof(Matrixf32MicrokernelTaskInfo) * number_of_rows);
		Matrixf32MicrokernelTaskInfo *this_task_info = microkernel_task_info_list;
		size_t task_number = 0;
		for (size_t row = 0; row < number_of_rows; row++) {
			// Each row defines a task. These are then run in parallel if there is a thread pool
			// otherwise they are run sequentially 
			/*
			auto rowTask = [=, &result]() {
				for (size_t col = 0; col < number_of_columns; col++) {
					matrixf32MicrokernelApply(result, a, b, row * kernel_height, col * kernel_width * SIMD_VECTOR_SIZE, kernel_width, kernel_height);
				}
				return true;
			};
			*/
			if (pool) {
				this_task_info->queue = &(pool->task_queue);
			}
			else {
				this_task_info->queue = nullptr;
			}
			this_task_info->kernel_width = kernel_width;
			this_task_info->kernel_height = kernel_height;
			this_task_info->left_matrix = a;
			this_task_info->right_matrix = b;
			this_task_info->result_row_number = row;
			this_task_info->result_matrix = result;
			this_task_info->task_number = task_number++;


			if (pool) {
				//TaskHandle rowFuture = pool->spawnTask(rowTask);
				//new(futures + row) TaskHandle(std::move(rowFuture));
				concurrentTaskQueuePush(this_task_info++, &matrixMicrokernelRowTask, &(pool->task_queue));
			}
			else {
				matrixMicrokernelRowTask(this_task_info++);
			}
		}
		for (size_t row = 0; row < remaining_rows_after_microkernel; row++) {

			auto rowTask = [=, &result]() {
				for (size_t col = 0; col < number_of_columns; col++) {
					matrixf32MicrokernelApply(result, a, b, number_of_rows * kernel_height + row, col * kernel_width * SIMD_VECTOR_SIZE, kernel_width, 1);
				}
				return true;
			};
			if (pool) {
				this_task_info->queue = &(pool->task_queue);
			}
			else {
				this_task_info->queue = nullptr;
			}
			this_task_info->kernel_width = kernel_width;
			this_task_info->kernel_height = 1;
			this_task_info->left_matrix = a;
			this_task_info->right_matrix = b;
			this_task_info->result_row_number = number_of_rows * kernel_height + row;
			this_task_info->result_matrix = result;
			this_task_info->task_number = task_number++;

			if (pool) {
				concurrentTaskQueuePush(this_task_info++, &matrixMicrokernelRowTask, &(pool->task_queue));
			//TaskHandle rowFuture = pool->spawnTask(rowTask);
			//new(futures + number_of_rows * kernel_height + row) TaskHandle(std::move(rowFuture));
			}
			else {
				matrixMicrokernelRowTask(this_task_info++);
			}
		}
		if (pool) {
			//pool->activeWait(*(futures + number_of_rows - 1));
			threadPoolActiveWait(pool);
			concurrentTaskQueueReset(&(pool->task_queue));
		}
		free(microkernel_task_info_list);
	}
	// Clean up what is left (columns)
	if (number_of_columns * kernel_width * SIMD_VECTOR_SIZE < b->cols)
	{
		for (size_t i = 0; i < a->rows; i++)
		{
			const float *aRow = a->contents + i * a->cols;
			float *result_row = result->contents + i * result->cols;
			for (size_t k = 0; k < b->rows; k++)
			{
				const float *bRow = b->contents + k * b->cols;
				const float Aik = aRow[k];
				for (size_t j = number_of_columns * kernel_width * SIMD_VECTOR_SIZE; j < b->cols; j++)
				{
					// is this compiler dependant? or just the overload issue?
#ifdef __APPLE__
					//todo(AION): same as above, this needs to be compiler dependant to be in line
					result_row[j] += Aik * bRow[j];
#else
					result_row[j] = fmaf(bRow[j], Aik, result_row[j]);
#endif

				}
			}
		}
	}
}

