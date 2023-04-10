#pragma once
#include <time.h>
#include "matrix.h"

#ifdef __APPLE__
typedef struct benchmark_endpoint {
	struct timespec time;
} benchmark_endpoint;

#else 
typedef struct benchmark_endpoint {
	clock_t benchmark_clock;
} benchmark_endpoint;
#endif
benchmark_endpoint GetBenchmarkEndpoint();
double GetTimeDiffFromEndpoints(benchmark_endpoint start, benchmark_endpoint end);

typedef struct benchmark_results
{
	double gflops;
	float max_error;
	size_t i;
	size_t j;
	size_t k;
	matrix_f32 result;
} benchmark_results;

void PrintBenchmark(benchmark_results benchmark);
double GFlops(double nOps, double timeInSeconds);
benchmark_results TestRandom(size_t size_i, size_t size_j, size_t size_k, int nTrials, const size_t kernelWidth, const size_t kernelHeight);
benchmark_results TestOnes(size_t size_i, size_t size_j, size_t size_k, int nTrials, const size_t kernelWidth, const size_t kernelHeight);
int MatrixBenchmarksMain(size_t sizeI, size_t sizeJ, size_t sizeK, int nTrials, const size_t kernelWidth, const size_t kernelHeight);
