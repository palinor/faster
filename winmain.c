#include <windows.h>
#include "tests.h"
#include "matrix_benchmarks.h"


/*int WINAPI wWinMain(
	HINSTANCE hInstance,
	HINSTANCE hPrevInstance,
	PWSTR lpCmdLine, 
	int nCmdShow) {
	TestsMain();
	printf("Hello world");
	return 0;
}*/

int main() {
	const size_t sizeI = 128;
	const size_t sizeJ = 128;
	const size_t sizeK = 128;
	const int nTrials = 1000;
	const size_t kernelWidth = 8;
	// todo(AION): Why does this segfault with kernelHeight = 3?
	const size_t kernelHeight = 2;
	MatrixBenchmarksMain(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight);
	return 0;
}