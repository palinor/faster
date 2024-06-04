#include <stdarg.h>
#include <stdlib.h>
#include "basic_math.h"
#include "linalg.h"
#include "matrix_benchmarks.h"
#include "date.cpp"

int LogErrorExact(float expectedOutput, float output, float *inputs, size_t nInputs)
{
	int error = !(ABS(expectedOutput - output) < 1e-8);
	if (error)
	{
		printf("Expected %lf, got %lf for inputs", expectedOutput, output);
		for (size_t i = 0; i < nInputs; i++)
		{
			printf(", %lf", inputs[i]);
		}
		printf("\n");
	}
	return error;
}

int LogErrorLeq(float expectedOutput, float output, float errorTol, float *inputs, size_t nInputs) {
	float absError = ABS(expectedOutput - output);
	int error = !(absError < errorTol);
	if (error)
	{
		printf("Expected diff of %lf and %lf to be under %lf, got %lf for inputs", output, expectedOutput, errorTol, absError);
		for (size_t i = 0; i < nInputs; i++) {
			printf(", %lf", inputs[i]);
		}
		printf("\n");
	}
	return error;
}

void RunTestFunction(int testF(void), const char *fName)

{
	int nErrors = testF();
	if (nErrors == 0)
	{
		printf("%s passed\n", fName);
	} else
	{
		printf("%s failed with %d errors\n", fName, nErrors);
	}
	printf("\n============\n");

}

int TestFSqrt()
{
	float inputs[] = {1, 0.001, 4, 100};
	float expectedOutputs[] = {1, 0.0316227749, 2.0, 10};
	printf("Testing FSqrt\n");
	int n_errors = 0;
	for (size_t i = 0; i < 4; i++)
	{
		float res = FSqrt(inputs[i]);
		n_errors += LogErrorExact(expectedOutputs[i], res, inputs + i, 1);
	}
	return n_errors;
}

int TestFPow()
{
	printf("Testing FPow\n");
	float baseInputs[] = {1, 10, 0.1, 3};
	int expInputs[] = {0, 1, 2, 5};
	float expectedRes[16] = {
		1, 1, 1, 1,
		1, 10, 100, 100000,
		1, 0.1, 0.01, 0.00001,
		1, 3, 9, 243
	};
	int nErrors = 0;
	for (size_t i = 0; i < 4; i++)
	{
		float base = baseInputs[i];
		for (size_t j = 0; j < 4; j++)
		{
			float exp = expInputs[j];
			float inputs[2] = {base, exp};
			float expected = expectedRes[4 * i + j];
			nErrors += LogErrorExact(expected, FPow(base, exp), inputs, 2);
		}
	}
	return nErrors;
}

void TestSwap() {
	float a = 5;
	float b = 6;
	Swap(&a, &b);
	printf("a: %f, b: %f", a, b);
}

/*
int TestGaussJordan() {
	printf("Testing GaussJordan \n");
	int nErrors = 0;
	size_t rowCases[4] = {8, 32, 128, 256};
	for (size_t testCase = 0; testCase < 4; testCase++) {
		size_t rows = rowCases[testCase];
		printf("Size %d\n", (int)rows);
		matrix_f32 id = Identityf32(rows);
		matrix_f32 a = RandomMatrixf32(rows, rows);
		matrix_f32 aInv = Matrixf32Copy(&a);
		GaussJordan(&aInv, NULL);
		matrix_f32 result;
		result.rows = a.rows;
		result.cols = aInv.cols;
		result.ld = result.cols;
		result.contents = (float *)calloc(result.rows * result.cols, sizeof(float));
		Matrixf32MicrokernelMultiply(&result, &a, &aInv, 4, 4, nullptr);
		float maxError = GetMaxError(&id, &result);
		float relErrorTol = 1e-6;
		float errorTol = ((float)rows) * relErrorTol;
		int doesNotPass = LogErrorLeq(0, maxError, errorTol, 0, 0);
		if (!doesNotPass) printf("Ok at %f\n", relErrorTol);
		nErrors += doesNotPass;
		if (rows < 10) {
			printf("Result: \n");
			PrintMatrixf32(&result);
		}
	}
	return nErrors;
}

int TestLU() {
	printf("Testing LU\n");
	int nErrors = 0;
	size_t rowCases[4] = {8, 32, 128, 256};
	for (size_t testCase = 0; testCase < 4; testCase++) {
		size_t rows = rowCases[testCase];
		printf("Size %d\n", (int)rows);
		matrix_f32 id = Identityf32(rows);
		matrix_f32 idCopy = Matrixf32Copy(&id);
		matrix_f32 a = RandomMatrixf32(rows, rows);
		matrix_f32 aLU = Matrixf32Copy(&a);
		LUFactorize(&aLU, 0);
		LUBackSub(&aLU, &id);
		// At this point our "id" matrix should be replaced with the inverse of A
		matrix_f32 result = Matrixf32Multiply(&a, &id);
		float relErrorTol = 1e-6;
		float errorTol = ((float)rows) * relErrorTol;
		float maxError = GetMaxError(&result, &idCopy);
		int doesNotPass = LogErrorLeq(0, maxError, errorTol, 0, 0);
		if (!doesNotPass) printf("Ok at %f\n", relErrorTol);
		nErrors += doesNotPass;
		if (rows < 10) {
			printf("Result: \n");
			PrintMatrixf32(&result);
		}
	}
	return nErrors;
}

void MatrixLUPopulate(matrix_f32 *lu, matrix_f32 *l, matrix_f32 *u) {
	for (size_t i = 0; i < l->rows; i++) {
		for (size_t j = 0; j < l->cols; j++) {
			if (i < j) {
				MATRIX_ITEM(u, i, j) = MATRIX_ITEM(lu, i, j);
				MATRIX_ITEM(l, i, j) = 0;
			}
			else if (i == j) {
				MATRIX_ITEM(u, i, j) = MATRIX_ITEM(lu, i, j);
				MATRIX_ITEM(l, i, j) = 1;
			}
			else {
				MATRIX_ITEM(u, i, j) = 0;
				MATRIX_ITEM(l, i, j) = MATRIX_ITEM(lu, i, j);
			}
		}
	}
}

matrix_f32 Matrixf32DumbMultiply(matrix_f32 *a, matrix_f32 *b) {
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = a->ld;
	result.contents = (float *)calloc(result.rows * result.cols, sizeof(float));
	for (size_t i = 0; i < a->rows; i++) {
		for (size_t j = 0; j < b->cols; j++) {
			float thisElem = 0;
			for (size_t k = 0; k < a->cols; k++) {
				thisElem += MATRIX_ITEM(a, i, k) * MATRIX_ITEM(b, k, j);
			}
			MATRIX_ITEM(&result, i, j) = thisElem;
		}
	}
	return result;
}

int TestLUFactorize() {
	size_t n = 6;
	int nErrors = 0;
	matrix_f32 l = LowerRandomf32(n);
	matrix_f32 u = UpperRandomf32(n);
	for (size_t i = 0; i < n; i++) {
		MATRIX_ITEM(&u, i, i) = 1;  // Enforce U to be unitary for solution uniqueness
	}
	matrix_f32 a;
	a.rows = n;
	a.cols = n;
	a.ld = n;
	a.contents = (float *)calloc(n * n, sizeof(float));
	// matrix_f32 a = Matrixf32Multiply(&l, &u);
	// Matrixf32MicrokernelMultiply(&a, &l, &u, 2, 2, nullptr);
	a = Matrixf32DumbMultiply(&l, &u);
	matrix_f32 aCopy = Matrixf32Copy(&a);
	size_t *index = reinterpret_cast<size_t*>(malloc(n * sizeof(size_t)));
	if (!index) {
		perror("Error allocating index in TestLUFactorize");
		exit(1);
	}
	LUFactorize(&a, index);
	MatrixLUPopulate(&a, &l, &u);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			Swap(MatrixGetAddr(&aCopy, index[i], j), MatrixGetAddr(&aCopy, i, j));
		}
	}
	matrix_f32 result;
	result.rows = n;
	result.cols = n;
	result.ld = n;
	result.contents = (float *)calloc(n * n, sizeof(float));
	// matrix_f32 result = Matrixf32Multiply(&l, &u);
	Matrixf32MicrokernelMultiply(&result, &l, &u, 2, 4, nullptr);
	float relErrorTol = 1e-8;
	float errorTol = ((float)n) * relErrorTol;
	float maxError = GetMaxError(&aCopy, &result);
	int doesNotPass = LogErrorLeq(0, maxError, errorTol, 0, 0);
		if (!doesNotPass) printf("Ok at %f\n", relErrorTol);
		nErrors += doesNotPass;
		if (n < 16) {
			printf("Result: \n");
			Matrixf32MultiplyInplace(&aCopy, -1);
			matrix_f32 error = Matrixf32Sum(&aCopy, &result);
			PrintMatrixf32(&error);
		}
	return nErrors;
}
*/

int TestCalendar() {
	Date start_date = {1, JANUARY, 2024};
	Date end_date = {2, JANUARY, 2024};
	printf("1 Jan 2024 - 2 Jan 2024, Calendar = None: %i", CoverageDays(&start_date, &end_date, Calendar::NONE));		
	printf("1 Jan 2024 - 2 Jan 2024, Calendar = Weekday: %i", CoverageDays(&start_date, &end_date, Calendar::WEEKDAYS));		
	printf("1 Jan 2024 - 2 Jan 2024, Calendar = SIFMA: %i", CoverageDays(&start_date, &end_date, Calendar::SIFMA));		
}

#ifdef IS_MAIN
/*
int main()
{
	const char fsqrtName[6] = "FSqrt";
	const char fpowName[5] = "FPow";
	const char gaussJordanName[12] = "GaussJordan";
	const char LUFactorizeName[13] = "LU Factorize";
	TestOnes(10, 10, 10, 1, 2, 4, nullptr);
	TestRandom(10, 10, 10, 1, 2, 4, nullptr);
	RunTestFunction(TestFSqrt, fsqrtName);
	RunTestFunction(TestFPow, fpowName);
	RunTestFunction(TestGaussJordan, gaussJordanName);
	RunTestFunction(TestLUFactorize, LUFactorizeName);
	// RunTestFunction(TestLU, "LUFactorize + LUBackSub");
	return 0;
} */
int main() {
	TestCalendar();
}
#endif
