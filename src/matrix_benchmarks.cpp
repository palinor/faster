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
	return TimeDiffInSeconds(start.time, end.time);
}

#else 
typedef struct benchmark_endpoint {
	std::chrono::high_resolution_clock::time_point benchmark_clock;;
} benchmark_endpoint;

benchmark_endpoint GetBenchmarkEndpoint() {
	benchmark_endpoint result;
	result.benchmark_clock = std::chrono::high_resolution_clock::now();
	return result;
}

double GetTimeDiffFromEndpoints(benchmark_endpoint start, benchmark_endpoint end) {
	std::chrono::duration<double> total_time = end.benchmark_clock - start.benchmark_clock;
	return std::chrono::duration_cast<std::chrono::nanoseconds>(total_time).count() / 1000000000.0;
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
	Matrixf32 result;
};

void PrintBenchmark(benchmark_results benchmark)
{
	printf("===========\n");
	printf("Matrix multiply %lu %lu %lu\n", benchmark.i, benchmark.j, benchmark.k);
	printf("GFLOPS: %lf\n", benchmark.gflops);
	printf("Max error: %lf\n", benchmark.max_error);
	printf("Total time: %lf\n", benchmark.total_time);
	printf("===========\n");
}


double GFlops(double nOps, double timeInSeconds)
{
	return nOps * 2e-9 / timeInSeconds;
}

float GetMaxError(Matrixf32 *a, Matrixf32 *b) {
	float result = 0;
	for (size_t i = 0; i < a->rows; ++i) {
		for (size_t j = 0; j < a->cols; ++j) {
			if (fabsf(matrixf32GetItem(a, i, j) - matrixf32GetItem(b, i, j)) > result) {
				result = fabsf(matrixf32GetItem(a, i, j) - matrixf32GetItem(b, i, j));
			}
		}
	}
	return result;
}

benchmark_results TestRandomNaive(
	size_t size_i, 
	size_t size_j, 
	size_t size_k, 
	int nTrials
)
{
	Matrixf32 a = matrixf32Random(size_i, size_k);
	Matrixf32 b = matrixf32Random(size_k, size_j);
	benchmark_results averageResult;
	averageResult.gflops = 0;
	averageResult.i = size_i;
	averageResult.j = size_j;
	averageResult.k = size_k;
	averageResult.total_time = 0;
	double nOps = size_i * size_j * size_k;
	Matrixf32 targetResult = matrixf32Multiply(&a, &b);
	for (int trial = 0; trial < nTrials; trial++)
	{
		benchmark_endpoint start = GetBenchmarkEndpoint();
		Matrixf32 testResult = matrixf32Multiply(&a, &b);
		benchmark_endpoint end = GetBenchmarkEndpoint();
		double timeInSeconds = GetTimeDiffFromEndpoints(start, end);
		averageResult.total_time += timeInSeconds;
		averageResult.max_error = GetMaxError(&testResult, &targetResult);
		averageResult.result = testResult;
		if (nTrials == 1) {
			matrixf32MultiplyInplace(&testResult, -1);
			Matrixf32 diff;
			matrixf32Sum(&diff, &targetResult, &testResult);
			matrixf32Print(&diff);
			matrixf32Free(&diff);
		}
		matrixf32Free(&testResult);
	}
	matrixf32Free(&a);
	matrixf32Free(&b);
	matrixf32Free(&targetResult);
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
	ThreadPool *pool
)
{
	Matrixf32 a = matrixf32Random(size_i, size_k);
	Matrixf32 b = matrixf32Random(size_k, size_j);
	benchmark_results averageResult;
	averageResult.gflops = 0;
	averageResult.i = size_i;
	averageResult.j = size_j;
	averageResult.k = size_k;
	averageResult.total_time = 0;
	double nOps = size_i * size_j * size_k;
	Matrixf32 targetResult = matrixf32Multiply(&a, &b);
	Matrixf32 testResult;
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
			matrixf32MultiplyInplace(&testResult, -1);
			Matrixf32 diff;
			matrixf32Sum(&diff, &targetResult, &testResult);
			matrixf32Print(&diff);
			matrixf32Free(&diff);
		}
	}
	averageResult.gflops = GFlops(nOps * nTrials, averageResult.total_time);
	matrixf32Free(&testResult);
	matrixf32Free(&targetResult);
	matrixf32Free(&a);
	matrixf32Free(&b);
	return averageResult;
}

benchmark_results TestOnes(
	size_t size_i, 
	size_t size_j, 
	size_t size_k, 
	int nTrials, 
	const size_t kernelWidth, 
	const size_t kernelHeight,
	ThreadPool *pool
)
{
	benchmark_results averageResult;
	Matrixf32 a = matrixf32Ones(size_i, size_k);
	Matrixf32 b = matrixf32Ones(size_k, size_j);
	averageResult.gflops = 0;
	averageResult.i = size_i;
	averageResult.j = size_j;
	averageResult.k = size_k;
	averageResult.total_time = 0;
	double nOps = size_i * size_j * size_k;
	Matrixf32 targetResult = matrixf32Multiply(&a, &b);
	Matrixf32 testResult;
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
			matrixf32MultiplyInplace(&testResult, -1);
			Matrixf32 diff;
			matrixf32Sum(&diff, &targetResult, &testResult);
			matrixf32Print(&diff);
		}
	}
	averageResult.gflops = GFlops(nOps * nTrials, averageResult.total_time);
	matrixf32Free(&a);
	matrixf32Free(&b);
	matrixf32Free(&testResult);
	matrixf32Free(&targetResult);
	return averageResult;
}

int MatrixBenchmarksMain(
	size_t sizeI, 
	size_t sizeJ, 
	size_t sizeK, 
	int nTrials, 
	const size_t kernelWidth, 
	const size_t kernelHeight,
	ThreadPool *pool
)
{

	printf("Testing Naive multiply on random matrices\n");
	benchmark_results benchmark = TestRandomNaive(sizeI, sizeJ, sizeK, nTrials);
	PrintBenchmark(benchmark);
	printf("Testing microkernel multiply on random matrices\n");
	benchmark = TestRandom(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight, pool);
	PrintBenchmark(benchmark);
	/*
	printf("Testing microkernel multiply on matrices of ones\n");
	benchmark = TestOnes(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight, pool);
	PrintBenchmark(benchmark);
	*/
	return 0;
}

int main() {
	const size_t sizeI = 256;
	const size_t sizeJ = 256;
	const size_t sizeK = 256;
	const int nTrials = 1000;
	const size_t kernelWidth = 4;
	const size_t kernelHeight = 4;
    //ThreadPool *pool = ThreadPool::getInstance();
	ThreadPool pool;
	//pool->start();
	size_t nStartingTasks = 256;
	//size_t nThreads = 15;
	threadPoolInit(&pool, nStartingTasks);
	printf("Testing C-style threadpool\n");
	MatrixBenchmarksMain(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight, &pool);
	//threadPoolInterrupt(&pool);
	return 0;
}
