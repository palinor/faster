#pragma once
#include <ctime>
#include <cstdlib>
#include <memory>
#include "matrix.cpp"

#ifdef __APPLE__
typedef struct benchmark_endpoint {
	struct timespec time;
} benchmark_endpoint;

struct timespec GetTime()
{
	struct timespec Time;
	clock_gettime(CLOCK_REALTIME, &Time);
	return Time;
}

benchmark_endpoint GetBenchmarkEndpoint() {
	benchmark_endpoint result;
	result.time = GetTime();
	return result;
}

double TimeDiffInSeconds(struct timespec start, struct timespec end)
{
	long seconds = end.tv_sec - start.tv_sec;
	long nanoseconds = end.tv_nsec - start.tv_nsec;
	return seconds + nanoseconds * 1e-9;
}

double GetTimeDiffFromEndpoints(benchmark_endpoint start, benchmark_endpoint end) {
	return TImeDiffInSeconds(start.time, end.time);
}

#else 
typedef struct benchmark_endpoint {
	clock_t benchmark_clock;
} benchmark_endpoint;

benchmark_endpoint GetBenchmarkEndpoint() {
	benchmark_endpoint result;
	result.benchmark_clock = clock();
	return result;
}

double GetTimeDiffFromEndpoints(benchmark_endpoint start, benchmark_endpoint end) {
	clock_t n_clocks = end.benchmark_clock - start.benchmark_clock;
	return ((double)n_clocks) / CLOCKS_PER_SEC;
}
#endif

struct benchmark_results
{
	double gflops;
	double total_time;
	float max_error;
	size_t i;
	size_t j;
	size_t k;
	matrix_f32 result;
};

void PrintBenchmark(benchmark_results benchmark)
{
	printf("===========\n");
	printf("Matrix multiply %llu %llu %llu\n", benchmark.i, benchmark.j, benchmark.k);
	printf("GFLOPS: %lf\n", benchmark.gflops);
	printf("Max error: %lf\n", benchmark.max_error);
	printf("Total time: %lf\n", benchmark.total_time);
	printf("===========\n");
}


double GFlops(double nOps, double timeInSeconds)
{
	return nOps * 2e-9 / timeInSeconds;
}

benchmark_results TestRandomNaive(
	size_t size_i, 
	size_t size_j, 
	size_t size_k, 
	int nTrials
)
{
	matrix_f32 a = RandomMatrixf32(size_i, size_k);
	matrix_f32 b = RandomMatrixf32(size_k, size_j);
	benchmark_results averageResult;
	averageResult.gflops = 0;
	averageResult.i = size_i;
	averageResult.j = size_j;
	averageResult.k = size_k;
	averageResult.total_time = 0;
	double nOps = size_i * size_j * size_k;
	matrix_f32 targetResult = Matrixf32Multiply(&a, &b);
	for (int trial = 0; trial < nTrials; trial++)
	{
		benchmark_endpoint start = GetBenchmarkEndpoint();
		matrix_f32 testResult = Matrixf32Multiply(&a, &b);
		benchmark_endpoint end = GetBenchmarkEndpoint();
		double timeInSeconds = GetTimeDiffFromEndpoints(start, end);
		averageResult.total_time += timeInSeconds;
		averageResult.max_error = GetMaxError(&testResult, &targetResult);
		averageResult.result = testResult;
		if (nTrials == 1) {
			Matrixf32MultiplyInplace(&testResult, -1);
			matrix_f32 diff = Matrixf32Sum(&targetResult, &testResult);
			PrintMatrixf32(&diff);
			FreeMatrixf32(&diff);
		}
		FreeMatrixf32(&testResult);
	}
	FreeMatrixf32(&a);
	FreeMatrixf32(&b);
	FreeMatrixf32(&targetResult);
	averageResult.gflops = GFlops(nOps * nTrials, averageResult.total_time);
	return averageResult;
}

benchmark_results TestRandom(
	size_t size_i, 
	size_t size_j, 
	size_t size_k, 
	int nTrials, 
	const size_t kernelWidth, 
	const size_t kernelHeight,
	thread_pool *pool
)
{
	matrix_f32 a = RandomMatrixf32(size_i, size_k);
	matrix_f32 b = RandomMatrixf32(size_k, size_j);
	benchmark_results averageResult;
	averageResult.gflops = 0;
	averageResult.i = size_i;
	averageResult.j = size_j;
	averageResult.k = size_k;
	averageResult.total_time = 0;
	double nOps = size_i * size_j * size_k;
	matrix_f32 targetResult = Matrixf32Multiply(&a, &b);
	matrix_f32 testResult;
	testResult.rows = a.rows;
	testResult.cols = b.cols;
	testResult.ld = b.ld;
	testResult.contents = (float *)calloc(a.rows * b.cols, sizeof(float));
	for (int trial = 0; trial < nTrials; trial++)
	{
		benchmark_endpoint start = GetBenchmarkEndpoint();
		Matrixf32MicrokernelMultiply(&testResult, &a, &b, kernelWidth, kernelHeight, pool);
		benchmark_endpoint end = GetBenchmarkEndpoint();
		double timeInSeconds = GetTimeDiffFromEndpoints(start, end);
		averageResult.total_time += timeInSeconds;
		averageResult.max_error = GetMaxError(&testResult, &targetResult);
		averageResult.result = testResult;
		averageResult.gflops += GFlops(nOps, timeInSeconds);
		if (nTrials == 1) {
			Matrixf32MultiplyInplace(&testResult, -1);
			matrix_f32 diff = Matrixf32Sum(&targetResult, &testResult);
			PrintMatrixf32(&diff);
			FreeMatrixf32(&diff);
		}
	}
	averageResult.gflops = GFlops(nOps * nTrials, averageResult.total_time);
	FreeMatrixf32(&testResult);
	FreeMatrixf32(&targetResult);
	FreeMatrixf32(&a);
	FreeMatrixf32(&b);
	return averageResult;
}

benchmark_results TestOnes(
	size_t size_i, 
	size_t size_j, 
	size_t size_k, 
	int nTrials, 
	const size_t kernelWidth, 
	const size_t kernelHeight,
	thread_pool *pool
)
{
	benchmark_results averageResult;
	matrix_f32 a = Onesf32(size_i, size_k);
	matrix_f32 b = Onesf32(size_k, size_j);
	averageResult.gflops = 0;
	averageResult.i = size_i;
	averageResult.j = size_j;
	averageResult.k = size_k;
	averageResult.total_time = 0;
	double nOps = size_i * size_j * size_k;
	matrix_f32 targetResult = Matrixf32Multiply(&a, &b);
	matrix_f32 testResult;
	testResult.rows = a.rows;
	testResult.cols = b.cols;
	testResult.ld = b.ld;
	testResult.contents = (float *)(calloc(a.rows * b.cols, sizeof(float)));
	for (int trial = 0; trial < nTrials; trial++)
	{
		benchmark_endpoint start = GetBenchmarkEndpoint();
		Matrixf32MicrokernelMultiply(&testResult, &a, &b, kernelWidth, kernelHeight, pool);
		benchmark_endpoint end = GetBenchmarkEndpoint();
		double timeInSeconds = GetTimeDiffFromEndpoints(start, end);
		averageResult.total_time += timeInSeconds;
		averageResult.max_error = GetMaxError(&testResult, &targetResult);
		averageResult.result = testResult;
		if (nTrials == 1) {
			Matrixf32MultiplyInplace(&testResult, -1);
			matrix_f32 diff = Matrixf32Sum(&targetResult, &testResult);
			PrintMatrixf32(&diff);
		}
	}
	averageResult.gflops = GFlops(nOps * nTrials, averageResult.total_time);
	FreeMatrixf32(&a);
	FreeMatrixf32(&b);
	FreeMatrixf32(&testResult);
	FreeMatrixf32(&targetResult);
	return averageResult;
}

int MatrixBenchmarksMain(
	size_t sizeI, 
	size_t sizeJ, 
	size_t sizeK, 
	int nTrials, 
	const size_t kernelWidth, 
	const size_t kernelHeight,
	thread_pool *pool
)
{
	//printf("Testing Naive multiply on random matrices\n");
	//benchmark_results benchmark = TestRandomNaive(sizeI, sizeJ, sizeK, nTrials);
	//PrintBenchmark(benchmark);
	printf("Testing microkernel multiply on random matrices\n");
	benchmark_results benchmark = TestRandom(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight, pool);
	PrintBenchmark(benchmark);
	printf("Testing microkernel multiply on matrices of ones\n");
	benchmark = TestOnes(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight, pool);
	PrintBenchmark(benchmark);
	return 0;
}
