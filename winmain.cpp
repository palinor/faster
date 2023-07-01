#include <windows.h>
#include "savinethreadpool.h"
#include "tests.h"
#include "matrix_benchmarks.h"

/*
int WINAPI wWinMain(
	HINSTANCE instance,
	HINSTANCE prevInstance,
	PWSTR commandLine,
	int showCode)
{
	WNDCLASS windowClass = {0};
	windowClass.style = CS_OWNDC | CS_HREDRAW | CS_VREDRAW;
	windowClass.hInstance = instance;
	windowClass.lpszClassName = "AnalyticsWindowClass";

	if (RegisterClass(&windowClass)) {
		HWND WindowHandle =
			CreateWindowEx();
	}
	return 0;
}
*/

int main() {
	const size_t sizeI = 256;
	const size_t sizeJ = 256;
	const size_t sizeK = 256;
	const int nTrials = 10000;
	const size_t kernelWidth = 4;
	const size_t kernelHeight = 4;
	ThreadPool *pool = ThreadPool::getInstance();
	pool->start(8);
	MatrixBenchmarksMain(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight);
	return 0;
}