#include <stdarg.h>
#include "basic_math.c"

int log_error(float expected_output, float output, float *inputs, size_t n_inputs)
{
	int error = !(ABS(expected_output - output) < 1e-8);
	if (error)
	{
		printf("Expected %lf, got %lf for inputs", expected_output, output);
		for (size_t i = 0; i < n_inputs; i++)
		{
			printf(", %lf", inputs[i]);
		}
		printf("\n");
	}
	return error;
}

void run_test_function(int test_f(void), char *f_name)

{
	int n_errors = test_f();
	if (n_errors == 0)
	{
		printf("%s passed\n", f_name);
	} else
	{
		printf("%s failed with %d errors\n", f_name, n_errors);
	}
	printf("\n============\n");

}

int test_f_sqrt()
{
	float inputs[] = {1, 0.001, 4, 100};
	float expected_outputs[] = {1, 0.0316227749, 2, 10};
	printf("Testing f_sqrt\n");
	int n_errors = 0;
	for (size_t i = 0; i < 4; i++)
	{
		float res = f_sqrt(inputs[i]);
		n_errors += log_error(expected_outputs[i], res, inputs + i, 1);
	}
	return n_errors;
}

int test_f_pow()
{
	printf("Testing f_pow\n");
	float base_inputs[] = {1, 10, 0.1, 3};
	int exp_inputs[] = {0, 1, 2, 5};
	float expected_res[16] = {
		1, 1, 1, 1,
		1, 10, 100, 100000,
		1, 0.1, 0.01, 0.00001,
		1, 3, 9, 243
	};
	int n_errors = 0;
	for (size_t i = 0; i < 4; i++)
	{
		float base = base_inputs[i];
		for (size_t j = 0; j < 4; j++)
		{
			float exp = exp_inputs[j];
			float inputs[2] = {base, exp};
			float expected = expected_res[4 * i + j];
			n_errors += log_error(expected, f_pow(base, exp), inputs, 2);
		}
	}
	return n_errors;
}

int main()
{
	run_test_function(test_f_sqrt, "f_sqrt");
	run_test_function(test_f_pow, "f_pow");
}
