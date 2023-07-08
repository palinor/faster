#include "matrix.cpp"

int main() {
	const size_t sizeI = 256;
	const size_t sizeJ = 256;
	const size_t sizeK = 256;
	const int nTrials = 1000;
	const size_t kernelWidth = 4;
	const size_t kernelHeight = 16;
	ThreadPool *pool = ThreadPool::getInstance();
	pool->start(10);
	MatrixBenchmarksMain(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight);
	return 0;
}
