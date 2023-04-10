#pragma once
#include <time.h>
#include "matrix.h"
#include "matrix_benchmarks.h"

#ifdef __APPLE__
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

void PrintBenchmark(benchmark_results benchmark)
{
	printf("===========\n");
	printf("Matrix multiply %llu %llu %llu\n", benchmark.i, benchmark.j, benchmark.k);
	printf("GFLOPS: %lf\n", benchmark.gflops);
	printf("Max error: %lf\n", benchmark.max_error);
	printf("===========\n");
}


double GFlops(double nOps, double timeInSeconds)
{
	return nOps * 2e-9 / timeInSeconds;
}


benchmark_results TestRandom(size_t size_i, size_t size_j, size_t size_k, int nTrials, const size_t kernelWidth, const size_t kernelHeight)
{
	const matrix_f32 a = RandomMatrixf32(size_i, size_k);
	const matrix_f32 b = RandomMatrixf32(size_k, size_j);
	benchmark_results averageResult;
	averageResult.gflops = 0;
	averageResult.i = a.rows;
	averageResult.j = b.cols;
	averageResult.k = a.cols;
	double nOps = averageResult.i * averageResult.j * averageResult.k;
	double time = 0;
	for (int trial = 0; trial < nTrials; trial++)
	{
		matrix_f32 targetResult = Matrixf32Multiply(&a, &b);
		benchmark_endpoint start = GetBenchmarkEndpoint();
		matrix_f32 testResult = Matrixf32MicrokernelMultiply(&a, &b, kernelWidth, kernelHeight);
		benchmark_endpoint end = GetBenchmarkEndpoint();
		float maxError = GetMaxError(&testResult, &targetResult);
		double timeInSeconds = GetTimeDiffFromEndpoints(start, end);
		time += timeInSeconds;
		averageResult.max_error = maxError;
		averageResult.result = testResult;
		if (nTrials == 1) {
			PrintMatrixf32(&targetResult);
			printf("\n\n\n");
			Matrixf32MultiplyInplace(&testResult, -1);
			matrix_f32 diff = Matrixf32Sum(&targetResult, &testResult);
			PrintMatrixf32(&diff);
		}
	}
	averageResult.gflops = nOps * (double)nTrials * 1e-9 / time;
	return averageResult;
}

benchmark_results TestOnes(size_t size_i, size_t size_j, size_t size_k, int nTrials, const size_t kernelWidth, const size_t kernelHeight)
{
	benchmark_results averageResult;
	const matrix_f32 a = Onesf32(size_i, size_k);
	const matrix_f32 b = Onesf32(size_k, size_j);
	averageResult.gflops = 0;
	averageResult.i = a.rows;
	averageResult.j = b.cols;
	averageResult.k = a.cols;
	for (int trial = 0; trial < nTrials; trial++)
	{
		benchmark_endpoint start = GetBenchmarkEndpoint();
		matrix_f32 testResult = Matrixf32MicrokernelMultiply(&a, &b, kernelWidth, kernelHeight);
		benchmark_endpoint end = GetBenchmarkEndpoint();
		matrix_f32 targetResult = Matrixf32Multiply(&a, &b);
		float maxError = GetMaxError(&testResult, &targetResult);
		double nOps = averageResult.i * averageResult.j * averageResult.k;
		double timeInSeconds = GetTimeDiffFromEndpoints(start, end);
		averageResult.gflops += GFlops(nOps, timeInSeconds);
		averageResult.max_error = maxError;
		averageResult.result = testResult;
	}
	averageResult.gflops /= (double)nTrials;

	return averageResult;
}

int MatrixBenchmarksMain(size_t sizeI, size_t sizeJ, size_t sizeK, int nTrials, const size_t kernelWidth, const size_t kernelHeight)
{
	benchmark_results benchmark = TestRandom(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight);
	PrintBenchmark(benchmark);
	benchmark = TestOnes(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight);
	PrintBenchmark(benchmark);
	return 0;
}
