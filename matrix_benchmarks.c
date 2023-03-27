#include "matrix.c"

struct timespec GetTime()
{
	struct timespec Time;
	clock_gettime(CLOCK_REALTIME, &Time);
	return Time;
}

double TimeDiffInSeconds(struct timespec start, struct timespec end)
{
	long seconds = end.tv_sec - start.tv_sec;
	long nanoseconds = end.tv_nsec - start.tv_nsec;
	return seconds + nanoseconds * 1e-9;
}

typedef struct benchmark_results
{
	double time;
	float max_error;
	size_t i;
	size_t j;
	size_t k;
	matrix_f32 result;
} benchmark_results;

double GFlops(benchmark_results benchmark)
{
	return benchmark.i * benchmark.j * benchmark.k * 2e-9 / benchmark.time;
}

void PrintBenchmark(benchmark_results benchmark)
{
	printf("===========\n");
	printf("Matrix multiply %li %li %li\n", benchmark.i, benchmark.j, benchmark.k);
	printf("GFLOPS: %lfs\n", GFlops(benchmark));
	printf("Max error: %lf\n", benchmark.max_error);
	printf("===========\n");
}


benchmark_results TestRandom(size_t size_i, size_t size_j, size_t size_k, matrix_f32 (*testFunction)(const matrix_f32 *, const matrix_f32 *), int nTrials)
{
	const matrix_f32 a = RandomMatrixf32(size_i, size_k);
	const matrix_f32 b = RandomMatrixf32(size_k, size_j);
	benchmark_results averageResult;
	averageResult.time = 0;
	averageResult.i = a.rows;
	averageResult.j = b.cols;
	averageResult.k = a.cols;
	for (int trial = 0; trial < nTrials; trial++)
	{
		matrix_f32 targetResult = Multiply(&a, &b);
		struct timespec startTime = GetTime();
		matrix_f32 testResult = testFunction(&a, &b);
		struct timespec endTime = GetTime();
		float maxError = GetMaxError(&testResult, &targetResult);
		averageResult.time += TimeDiffInSeconds(startTime, endTime);
		averageResult.max_error = maxError;
		averageResult.result = testResult;
	}
	averageResult.time /= (double)nTrials;

	return averageResult;
}

benchmark_results TestOnes(size_t size_i, size_t size_j, size_t size_k, matrix_f32 (*testFunction)(const matrix_f32 *, const matrix_f32 *), int nTrials)
{
	benchmark_results averageResult;
	const matrix_f32 a = Onesf32(size_i, size_k);
	const matrix_f32 b = Onesf32(size_k, size_j);
	averageResult.time = 0;
	averageResult.i = a.rows;
	averageResult.j = b.cols;
	averageResult.k = a.cols;
	for (int trial = 0; trial < nTrials; trial++)
	{
		struct timespec startTime = GetTime();
		matrix_f32 testResult = testFunction(&a, &b);
		struct timespec endTime = GetTime();
		matrix_f32 targetResult = Multiply(&a, &b);
		float maxError = GetMaxError(&testResult, &targetResult);
		averageResult.time += TimeDiffInSeconds(startTime, endTime);
		averageResult.max_error = maxError;
		averageResult.result = testResult;
	}
	averageResult.time /= (double)nTrials;

	return averageResult;
}

int main()
{
	size_t size_i = 128;
	size_t size_j = 128;
	size_t size_k = 128;
	int nTrials = 1000;
	benchmark_results benchmark = TestRandom(size_i, size_j, size_k, &MatrixMicrokernelMultiply, nTrials);
	PrintBenchmark(benchmark);
	benchmark = TestOnes(size_i, size_j, size_k, &MatrixMicrokernelMultiply, nTrials);
	PrintBenchmark(benchmark);
}
