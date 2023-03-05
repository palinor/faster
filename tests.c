#include <stdarg.h>
#include "basic_math.c"
#include "linalg.c"

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

void RunTestFunction(int testF(void), char *fName)

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
	float expectedOutputs[] = {1, 0.0316227749, 2, 10};
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
	swap(&a, &b);
	printf("a: %f, b: %f", a, b);
}

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
		int error = GaussJordan(&aInv, NULL);
		matrix_f32 result = Matrixf32Multiply(&a, &aInv);
		float maxError = GetMaxError(&id, &result);
		float relErrorTol = 1e-6;
		float errorTol = ((float)rows) * relErrorTol;
		int doesNotPass = LogErrorLeq(0, maxError, errorTol, 0, 0);
		if (!doesNotPass) printf("Ok at %f\n", relErrorTol);
		nErrors += doesNotPass;
		if (rows < 10) {
			printf("Result: \n");
			PrintMatrix(&result);
		}
	}
	return nErrors;
}

int main()
{
	RunTestFunction(TestFSqrt, "FSqrt");
	RunTestFunction(TestFPow, "FPow");
	RunTestFunction(TestGaussJordan, "GaussJordan");

}
