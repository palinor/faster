#include <windows.h>
//#include "matrix_benchmarks.cpp"
#include "aad.cpp"
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
/*
int main() {
	const size_t sizeI = 256;
	const size_t sizeJ = 256;
	const size_t sizeK = 256;
	const int nTrials = 100000;
	const size_t kernelWidth = 4;
	const size_t kernelHeight = 4;
	thread_pool pool;
	size_t nStartingTasks = 256;
	size_t nThreads = 14;
	InitThreadPool(&pool, nStartingTasks, nThreads);
	MatrixBenchmarksMain(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight, &pool);
	StopThreadPool(&pool);
	return 0;
}
*/
int main() {
	TestShiftedSabr();
}