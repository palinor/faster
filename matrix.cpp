#ifdef __APPLE__
#include <arm_neon.h>
#else 
#include <immintrin.h>
#include <omp.h>
#endif

#ifdef _WIN32
#pragma warning( disable: 4305 4244 6386)
#endif

#include <cassert>
#include <cmath>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <vector>

#include <chrono>
#include <condition_variable>
#include <future>
#include <queue>
#include <mutex>
#include <thread>
#include <stdio.h>

#define ABS(x) (((x) > 0) ? (x) : -(x))
#define MAX_FLOAT 1e20
#define MIN_FLOAT 1e-12
#define MATRIX_ITEM(a, i, j) ((a)->contents[(i) * (a)->ld + (j)])
#define SIGN(x) (((x) > 0) ? 1 : -1)


using task = std::packaged_task<bool(void)>;
using task_handle = std::future<bool>;

struct task_queue {
	task *tasks;
	size_t startingMaxSize;
	size_t maxSize;
	size_t currentSize;
	size_t currentHead;
};

void CreateTaskQueue(task_queue *result, size_t maxSize) {
	result->tasks = (task *)malloc(maxSize * sizeof(task));
	result->maxSize = maxSize;
	result->startingMaxSize = maxSize;
	result->currentSize = 0;
	result->currentHead = 0;
}

void FreeTaskQueue(task_queue *result) {
	free(result->tasks);
}

void TaskQueuePush(task_queue *input, task *element) {
	if (input->currentSize == input->maxSize) {
		size_t newSize = input->maxSize * 2 + 1;
		input->tasks = (task *)realloc(input->tasks, newSize * sizeof(task));
		input->maxSize = newSize;
	}
	input->tasks[input->currentSize++] = std::move(*element);
}

inline bool IsTaskQueueEmpty(task_queue *input) {
	return (input->currentHead == input->maxSize);
}

int TaskQueuePop(task *output, task_queue *input) {
	if (input->currentHead == input->maxSize) {
		return 1;
	}
	*output = std::move(input->tasks[input->currentHead]);
	input->currentHead++;

	if (IsTaskQueueEmpty(input)) {
		size_t newSize = input->maxSize - input->startingMaxSize;
		task *oldTasks = input->tasks;
		task *newTasks = (task *)realloc(oldTasks + input->startingMaxSize, newSize * sizeof(task));
		input->tasks = newTasks;
		free(oldTasks);
	}
	return 0;
}


struct concurrent_task_queue {
	task_queue queue_;
	mutable std::mutex mutex_;
	std::condition_variable cv_;
	bool interrupt_;
};

void InterruptConcurrentTaskQueue(concurrent_task_queue *input) {
	input->mutex_.lock();
	input->interrupt_ = true;
	input->mutex_.unlock();
	input->cv_.notify_all();
}

inline bool IsConcurrentTaskQueueEmpty(concurrent_task_queue *input) {
	return IsTaskQueueEmpty(&(input->queue_));
}

void ConcurrentTaskQueuePush(concurrent_task_queue *input, task *newTask) {
	input->mutex_.lock();
	TaskQueuePush(&(input->queue_), newTask);
	input->mutex_.unlock();
	input->cv_.notify_one();
}

bool ConcurrentTaskQueuePop(task *result, concurrent_task_queue *input) {
	std::unique_lock<std::mutex> lk(input->mutex_);
	while (!input->interrupt_ && IsTaskQueueEmpty(&(input->queue_))) {
		input->cv_.wait(lk);
	}
	if (input->interrupt_) {
		return false;
	}
	TaskQueuePop(result, &(input->queue_));
	return true;
}

bool ConcurrentTaskQueuePopTry(task *result, concurrent_task_queue *input) {
	input->mutex_.lock();
	if (IsTaskQueueEmpty(&(input->queue_))) {
		return false;
	}
	TaskQueuePop(result, &(input->queue_));
	input->mutex_.unlock();
	return true;
}

bool tryPop(T &t) {
	std::lock_guard<std::mutex> lk(my_mutex_);
	if (my_queue_.empty()) return false;
	t = std::move(my_queue_.front());
	my_queue_.pop();
	return true;
}

void clear() {
	std::queue<T> empty;
	swap(my_queue_, empty);
}

template <class T>
class ConcurrentQueue {
	std::queue<T> my_queue_;
	mutable std::mutex my_mutex_;
	std::condition_variable my_cv_;
	bool my_interrupt_;

public:

	ConcurrentQueue() : my_interrupt_(false) {}
	~ConcurrentQueue() { interrupt(); }

	void interrupt() {
		{
			std::lock_guard<std::mutex> lk(my_mutex_);
			my_interrupt_ = true;
		}
		my_cv_.notify_all();
	}

	void resetInterrupt() {
		my_interrupt_ = false;
	}

	bool empty() const {
		std::lock_guard<std::mutex> lk(my_mutex_);
		return my_queue_.empty();
	}

	void push(T t) {
		std::lock_guard<std::mutex> lk(my_mutex_);
		my_queue_.push(std::move(t));
		my_cv_.notify_one();
	}

	bool pop(T &t) {
		std::unique_lock<std::mutex> lk(my_mutex_);
		while (!my_interrupt_ && my_queue_.empty()) my_cv_.wait(lk);
		if (my_interrupt_) return false;
		t = std::move(my_queue_.front());
		my_queue_.pop();
		return true;
	}

	bool tryPop(T &t) {
		std::lock_guard<std::mutex> lk(my_mutex_);
		if (my_queue_.empty()) return false;
		t = std::move(my_queue_.front());
		my_queue_.pop();
		return true;
	}

	void clear() {
		std::queue<T> empty;
		swap(my_queue_, empty);
	}
};


typedef struct ThreadPool {
	static thread_local size_t my_thread_number_;


} ThreadPool;

void RunSingleThread(ConcurrentQueue<Task> *workQueue, bool *interrupt) {
	Task t;
	while (!interrupt) {
		workQueue->pop(t);
		if (!interrupt) t();
	}
}

class ThreadPool {

	static ThreadPool my_instance_;
	static thread_local size_t my_thread_number_;
	std::vector<std::thread> my_threads_;
	bool is_active_;
	bool is_interrupt_;
	ConcurrentQueue<Task> my_queue_;
	ThreadPool() : is_active_(false), is_interrupt_(false) {}

	void threadFunc(const size_t thread_number);

public:
	static ThreadPool *getInstance() { return &my_instance_; }

	ThreadPool(const ThreadPool &rhs) = delete;
	ThreadPool &operator=(const ThreadPool &rhs) = delete;
	ThreadPool(ThreadPool &&rhs) = delete;
	ThreadPool &operator=(ThreadPool &&rhs) = delete;

	void start(const size_t n_thread = std::thread::hardware_concurrency() - 1);
	size_t numThreads() const { return my_threads_.size(); }
	static size_t threadNum() { return my_thread_number_; }

	void stop();
	~ThreadPool() {
		stop();
	}

	template<typename Callable>
	TaskHandle spawnTask(Callable c) {
		Task t(std::move(c));
		TaskHandle f = t.get_future();
		my_queue_.push(std::move(t));
		return f;
	}

	bool activeWait(const TaskHandle &f);

};


using namespace std::chrono_literals;

void ThreadPool::threadFunc(const size_t thread_number) {
	my_thread_number_ = thread_number;
	Task t;
	while (!is_interrupt_) {
		my_queue_.pop(t);
		if (!is_interrupt_) t();
	}
}

void ThreadPool::start(const size_t n_thread) {
	if (!is_active_) {
		my_threads_.reserve(n_thread);
		for (size_t i = 0; i < n_thread; i++)
			my_threads_.push_back(std::thread(&ThreadPool::threadFunc, this, i + 1));
		is_active_ = true;
	}
}


void ThreadPool::stop() {
	if (is_active_) {
		is_interrupt_ = true;
		my_queue_.interrupt();
		std::for_each(
			my_threads_.begin(),
			my_threads_.end(),
			std::mem_fn(&std::thread::join)
		);
		my_threads_.clear();
		my_queue_.clear();
		my_queue_.resetInterrupt();
		is_active_ = false;
		is_interrupt_ = false;
	}
}


bool ThreadPool::activeWait(const TaskHandle &f) {
	Task t;
	bool i_did_work = false;
	while (f.wait_for(0s) != std::future_status::ready) {
		if (my_queue_.tryPop(t)) {
			t();
			i_did_work = true;
		} else {
			f.wait();
		}
	}
	return i_did_work;
}


ThreadPool ThreadPool::my_instance_;
thread_local size_t ThreadPool::my_thread_number_ = 0;

typedef struct matrix_f32
{
	size_t rows;
	size_t cols;
	float *contents;
	size_t ld;
} matrix_f32;

void FreeMatrixf32(matrix_f32 *a) {
	free(a->contents);
}

matrix_f32 Matrixf32Copy(matrix_f32 *a) {
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = a->cols;
	result.ld = a->ld;
	result.contents = reinterpret_cast<float *>(malloc(a->rows * a->cols * sizeof(float)));
	if (!result.contents) {
		perror("Error allocating result.contents in Matrixf32Copy");
		exit(1);
	}
	for (size_t i = 0; i < a->rows; i++) {
		for (size_t j = 0; j < a->cols; j++) {
			MATRIX_ITEM(&result, i, j) = MATRIX_ITEM(a, i, j);
		}
	}
	return result;
}

float MatrixGetItem(matrix_f32 *a, size_t i, size_t j)
{
	return a->contents[i * a->ld + j];
}

void MatrixPutItem(matrix_f32 *a, size_t i, size_t j, float x) {
	a->contents[i * a->ld + j] = x;
}

float *MatrixGetAddr(const matrix_f32 *a, size_t i, size_t j)
{
	return a->contents + i * a->ld + j;
}


void PrintMatrixf32(matrix_f32 *result)
{
	for (size_t i = 0; i < result->rows; i++)
	{
		for (size_t j = 0; j < result->cols; j++)
		{
			printf("%f ", MatrixGetItem(result, i, j));
		}
		printf("\n");
	}
}

double GetMaxError(matrix_f32 *a, matrix_f32 *b)
{
	assert(a->rows == b->rows);
	assert(a->cols == b->cols);
	float maxError = 0;
	for (size_t i = 0; i < a->rows; i++)
	{
		for (size_t j = 0; j < a->cols; j++)
		{
			float thisError = ABS(MatrixGetItem(a, i, j) - MatrixGetItem(b, i, j));
			if (!(thisError < MAX_FLOAT)) {
				return MAX_FLOAT; // handle NaN
			}
			maxError = (thisError > maxError) ? thisError : maxError;
		}
	}
	return maxError;
}


int Matrixf32MultiplyToTarget(matrix_f32 *result, const matrix_f32 *a, const matrix_f32 *b)
{
	assert(a->cols == b->rows);
	// We check that result has already been allocated properly
	assert(result->rows == a->rows);
	assert(result->cols == b->cols);
	assert(result->ld == b->cols);
	assert(!!(result->contents));
	if (!result->contents) {
		perror("Error allocating result.contents in Matrixf32Multiply");
		exit(1);
	}
	for (size_t i = 0; i < a->rows; i++)
	{
		const float *aRow = a->contents + i * a->cols;
		float *resultRow = result->contents + i * result->cols;
		for (size_t k = 0; k < b->rows; k++)
		{
			const float *bRow = b->contents + k * b->cols;
			const float Aik = aRow[k];
			// #pragma clang loop vectorize_width(8) // interleave_count(4)
			for (size_t j = 0; j < b->cols; j++)
			{
				//todo(AION): Have this swap depending on the platform we compile on.
				//resultRow[j] += Aik * bRow[j];
				/*
				This part is actually out of line with the microkernel loop (no difference in -Ofast, but visible in other compile modes)
				The above lines up just fine in all

				BUUUUUUUT.... not on windows. So we actually need to swap depending on the compiler. Lol.
				*/
				resultRow[j] = fma(bRow[j], Aik, resultRow[j]);
			}
		}
	}
	return 0;
}

matrix_f32 Matrixf32Multiply(const matrix_f32 *a, const matrix_f32 *b)
{
	assert(a->cols == b->rows);
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = b->cols;
	result.contents = reinterpret_cast<float *>(calloc(result.rows * result.cols, sizeof(float)));
	int err = Matrixf32MultiplyToTarget(&result, a, b);
	return result;
}
/*
Multiply inplace by a scalar
*/
void Matrixf32MultiplyInplace(matrix_f32 *a, float x) {
	for (size_t i = 0; i < a->rows; i++) {
		for (size_t j = 0; j < a->cols; j++) {
			MATRIX_ITEM(a, i, j) *= x;
		}
	}
}

matrix_f32 Matrixf32Sum(matrix_f32 *a, matrix_f32 *b) {
	assert(a->rows == b->rows);
	assert(a->cols == b->cols);
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = result.cols;
	result.contents = reinterpret_cast<float *>(malloc(a->rows * a->cols * sizeof(float)));
	if (!result.contents) {
		perror("Error allocating results in Matrixf32Sum");
		exit(1);
	}
	for (size_t i = 0; i < result.rows; i++) {
		for (size_t j = 0; j < result.cols; j++) {
			result.contents[i * result.ld + j] = MATRIX_ITEM(a, i, j) + MATRIX_ITEM(b, i, j);
		}
	}
	return result;
}


matrix_f32 RandomMatrixf32(size_t rows, size_t cols)
{
	matrix_f32 result;
	result.rows = rows;
	result.cols = cols;
	result.ld = cols;
	result.contents = reinterpret_cast<float *>(malloc(rows * cols * sizeof(float)));
	if (!result.contents) {
		perror("Error allocating result.contents in RandomMatrixf32");
		exit(1);
	}
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			result.contents[i * cols + j] = ((float)rand()) / RAND_MAX;
		}
	}
	return result;
}

matrix_f32 Onesf32(size_t rows, size_t cols)
{
	matrix_f32 result;
	result.rows = rows;
	result.cols = cols;
	result.ld = cols;
	result.contents = reinterpret_cast<float *>(malloc(rows * cols * sizeof(float)));
	if (!result.contents) {
		perror("Error allocating result.contents in Onesf32");
		exit(1);
	}

	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			result.contents[i * cols + j] = 1;
		}
	}
	return result;
}

matrix_f32 Identityf32(size_t n)
{
	matrix_f32 result;
	result.rows = n;
	result.cols = n;
	result.ld = n;
	result.contents = reinterpret_cast<float *>(calloc(n * n, sizeof(float)));
	for (size_t i = 0; i < n; i++)
	{
		result.contents[i * result.ld + i] = 1;
	}
	return result;
}

/*
Returns an upper triangular matrix with random non-zero elements
*/
matrix_f32 UpperRandomf32(size_t n) {
	matrix_f32 result;
	result.rows = n;
	result.cols = n;
	result.ld = n;
	result.contents = reinterpret_cast<float *>(calloc(n * n, sizeof(float)));
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i; j < n; j++) {
			result.contents[i * result.ld + j] = 100 * ((float)rand()) / RAND_MAX;
		}
	}
	return result;
}

/*
Returns an lower triangular matrix with random non-zero elements
*/
matrix_f32 LowerRandomf32(size_t n) {
	matrix_f32 result;
	result.rows = n;
	result.cols = n;
	result.ld = n;
	result.contents = reinterpret_cast<float *>(calloc(n * n, sizeof(float)));
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j > 0; j--) { // stupid trick because j is an unsigned int
			MATRIX_ITEM(&result, i, j - 1) = 100 * ((float)rand()) / RAND_MAX;
		}
	}
	return result;
}

/*
Internal block for matrix multiplication - compute the result in blocks of isze kernelHeight x (4 * kernelWidth)
*/
#ifdef __APPLE__
// ARM implimentation
#define SIMD_VECTOR_SIZE 4
void Microkernel(const matrix_f32 *a, const matrix_f32 *b, float *resultLocation, size_t resultRow, size_t resultCol, size_t kernelWidth, size_t kernelHeight)
{
	assert(a->cols == b->rows);
	assert(kernelWidth * kernelHeight < 128);
	// resultElems is really an array of size kernelWidth * kernelHeight
	float32x4_t resultElems[128];
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			resultElems[i * kernelWidth + j] = vdupq_n_f32(0);
		}
	}
	for (size_t k = 0; k < a->cols / SIMD_VECTOR_SIZE; k++)
	{
		for (size_t i = 0; i < kernelHeight; i++)
		{
			float32x4_t aElems_i = vld1q_f32(MatrixGetAddr(a, resultRow + i, SIMD_VECTOR_SIZE * k));
			for (size_t j = 0; j < kernelWidth; j++)
			{
				float32x4_t bElems_i0 = vld1q_f32(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k, resultCol + SIMD_VECTOR_SIZE * j));
				resultElems[i * kernelWidth + j] = vmlaq_laneq_f32(resultElems[i * kernelWidth + j], bElems_i0, aElems_i, 0);

				float32x4_t bElems_i1 = vld1q_f32(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 1, resultCol + SIMD_VECTOR_SIZE * j));
				resultElems[i * kernelWidth + j] = vmlaq_laneq_f32(resultElems[i * kernelWidth + j], bElems_i1, aElems_i, 1);

				float32x4_t bElems_i2 = vld1q_f32(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 2, resultCol + SIMD_VECTOR_SIZE * j));
				resultElems[i * kernelWidth + j] = vmlaq_laneq_f32(resultElems[i * kernelWidth + j], bElems_i2, aElems_i, 2);

				float32x4_t bElems_i3 = vld1q_f32(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 3, resultCol + SIMD_VECTOR_SIZE * j));
				resultElems[i * kernelWidth + j] = vmlaq_laneq_f32(resultElems[i * kernelWidth + j], bElems_i3, aElems_i, 3);
			}
		}
	}

	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			vst1q_f32(resultLocation + i * b->cols + j * SIMD_VECTOR_SIZE, resultElems[i * kernelWidth + j]);
		}
	}
}

/*
If k is not divisible by SIMD_VECTOR_SIZE, handle the addition of the remaining elements
*/

void MicrokernelRemainder(const matrix_f32 *a, const matrix_f32 *b, float *resultLocation, size_t resultRow, size_t resultCol, size_t kernelWidth, size_t kernelHeight)
{
	size_t remainder = a->cols % SIMD_VECTOR_SIZE;
	float32x4_t resultElems[32];
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			resultElems[i * kernelWidth + j] = vdupq_n_f32(0);
		}
	}
	for (size_t i = 0; i < kernelHeight; i++)
	{
		float32x4_t aRemainder_i = vld1q_f32(MatrixGetAddr(a, resultRow + i, a->cols - remainder));
		for (size_t j = 0; j < kernelWidth; j++)
		{
			float32x4_t bRemainderElems_0j = vld1q_f32(MatrixGetAddr(b, a->cols - remainder, resultCol + SIMD_VECTOR_SIZE * j));
			resultElems[i * kernelWidth + j] = vmlaq_laneq_f32(resultElems[i * kernelWidth + j], bRemainderElems_0j, aRemainder_i, 0);
			if (remainder > 1)
			{
				float32x4_t bRemainderElems_1j = vld1q_f32(MatrixGetAddr(b, a->cols - remainder + 1, resultCol + SIMD_VECTOR_SIZE * j));
				resultElems[i * kernelWidth + j] = vmlaq_laneq_f32(resultElems[i * kernelWidth + j], bRemainderElems_1j, aRemainder_i, 1);
			}
			if (remainder > 2)
			{
				float32x4_t bRemainderElems_2j = vld1q_f32(MatrixGetAddr(b, a->cols - remainder + 2, resultCol + SIMD_VECTOR_SIZE * j));
				resultElems[i * kernelWidth + j] = vmlaq_laneq_f32(resultElems[i * kernelWidth + j], bRemainderElems_2j, aRemainder_i, 2);
			}
		}
	}
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			vst1q_f32(resultLocation + i * b->cols + j * SIMD_VECTOR_SIZE, resultElems[i * kernelWidth + j]);
		}
	}
}
#else
// AMD64 implementation
#define SIMD_VECTOR_SIZE 8 
void Microkernel(const matrix_f32 *a, const matrix_f32 *b, float *resultLocation, size_t resultRow, size_t resultCol, const size_t kernelWidth, const size_t kernelHeight)
{
	assert(a->cols == b->rows);
	assert(kernelHeight * kernelWidth < 32);
	__m256 resultElems[32];
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			resultElems[i * kernelWidth + j] = _mm256_set1_ps(0);
		}
	}
	for (size_t k = 0; k < a->cols / SIMD_VECTOR_SIZE; k++)
	{
		for (size_t i = 0; i < kernelHeight; i++)
		{
			float *aElems_i = MatrixGetAddr(a, resultRow + i, SIMD_VECTOR_SIZE * k);
			for (size_t j = 0; j < kernelWidth; j++)
			{
				//todo(AION): we have 3 fma ports, why can't we use 3 accumulators? It seemed needlessly slow when I tried 2
				__m256 bElems_0j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i0 = _mm256_set1_ps(*aElems_i);
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i0, bElems_0j, resultElems[i * kernelWidth + j]);

				__m256 bElems_1j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 1, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i1 = _mm256_set1_ps(*(aElems_i + 1));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i1, bElems_1j, resultElems[i * kernelWidth + j]);

				__m256 bElems_2j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 2, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i2 = _mm256_set1_ps(*(aElems_i + 2));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i2, bElems_2j, resultElems[i * kernelWidth + j]);

				__m256 bElems_3j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 3, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i3 = _mm256_set1_ps(*(aElems_i + 3));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i3, bElems_3j, resultElems[i * kernelWidth + j]);

				__m256 bElems_4j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 4, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i4 = _mm256_set1_ps(*(aElems_i + 4));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i4, bElems_4j, resultElems[i * kernelWidth + j]);

				__m256 bElems_5j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 5, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i5 = _mm256_set1_ps(*(aElems_i + 5));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i5, bElems_5j, resultElems[i * kernelWidth + j]);

				__m256 bElems_6j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 6, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i6 = _mm256_set1_ps(*(aElems_i + 6));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i6, bElems_6j, resultElems[i * kernelWidth + j]);

				__m256 bElems_7j = _mm256_load_ps(MatrixGetAddr(b, SIMD_VECTOR_SIZE * k + 7, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 broadcastA_i7 = _mm256_set1_ps(*(aElems_i + 7));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(broadcastA_i7, bElems_7j, resultElems[i * kernelWidth + j]);
			}
		}
	}

	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			_mm256_store_ps(resultLocation + i * b->cols + j * SIMD_VECTOR_SIZE, resultElems[i * kernelWidth + j]);
		}
	}

}

/*
If k is not divisible by SIMD_VECTOR_SIZE, handle the addition of the remaining elements
*/

void MicrokernelRemainder(const matrix_f32 *a, const matrix_f32 *b, matrix_f32 *result, size_t resultRow, size_t resultCol, size_t kernelWidth, size_t kernelHeight)
{
	size_t remainder = a->cols % SIMD_VECTOR_SIZE;
	__m256 resultElems[32]:
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			resultElems[i * kernelWidth + j] = _mm256_load_ps(MatrixGetAddr(result, resultRow + i, resultCol + SIMD_VECTOR_SIZE * j));
		}
	}
	for (size_t i = 0; i < kernelHeight; i++)
	{
		float *aRemainder_i = MatrixGetAddr(a, resultRow + i, a->cols - remainder);
		for (size_t j = 0; j < kernelWidth; j++)
		{
			__m256 bRemainderElems_0j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder, resultCol + SIMD_VECTOR_SIZE * j));
			__m256 aRemainderElems_i0 = _mm256_set1_ps(*aRemainder_i);
			resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i0, bRemainderElems_0j, resultElems[i * kernelWidth + j]);
			if (remainder > 1)
			{
				__m256 bRemainderElems_1j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 1, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i1 = _mm256_set1_ps(*(aRemainder_i + 1));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i1, bRemainderElems_1j, resultElems[i * kernelWidth + j]);
			}
			if (remainder > 2)
			{
				__m256 bRemainderElems_2j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 2, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i2 = _mm256_set1_ps(*(aRemainder_i + 2));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i2, bRemainderElems_2j, resultElems[i * kernelWidth + j]);
			}
			if (remainder > 3)
			{
				__m256 bRemainderElems_3j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 3, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i3 = _mm256_set1_ps(*(aRemainder_i + 3));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i3, bRemainderElems_3j, resultElems[i * kernelWidth + j]);
			}
			if (remainder > 4)
			{
				__m256 bRemainderElems_4j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 4, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i4 = _mm256_set1_ps(*(aRemainder_i + 4));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i4, bRemainderElems_4j, resultElems[i * kernelWidth + j]);
			}
			if (remainder > 5)
			{
				__m256 bRemainderElems_5j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 5, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i5 = _mm256_set1_ps(*(aRemainder_i + 5));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i5, bRemainderElems_5j, resultElems[i * kernelWidth + j]);
			}
			if (remainder > 6)
			{
				__m256 bRemainderElems_6j = _mm256_load_ps(MatrixGetAddr(b, a->cols - remainder + 6, resultCol + SIMD_VECTOR_SIZE * j));
				__m256 aRemainderElems_i6 = _mm256_set1_ps(*(aRemainder_i + 6));
				resultElems[i * kernelWidth + j] = _mm256_fmadd_ps(aRemainderElems_i6, bRemainderElems_6j, resultElems[i * kernelWidth + j]);
			}
		}
	}
	for (size_t i = 0; i < kernelHeight; i++)
	{
		for (size_t j = 0; j < kernelWidth; j++)
		{
			_mm256_store_ps(MatrixGetAddr(result, resultRow + i, resultCol + SIMD_VECTOR_SIZE * j), resultElems[i * kernelWidth + j]);
		}
	}
}
#endif
matrix_f32 Matrixf32MicrokernelMultiply(const matrix_f32 *a, const matrix_f32 *b, const size_t kernelWidth, const size_t kernelHeight)
{
	assert(a->cols == b->rows);
	matrix_f32 result;
	result.rows = a->rows;
	result.cols = b->cols;
	result.ld = b->cols;
	result.contents = reinterpret_cast<float *>(calloc(a->rows * b->cols, sizeof(float)));
	const size_t nRows = a->rows / kernelHeight;
	const size_t rowRemainder = a->rows % kernelHeight;
	const size_t nCols = b->cols / (kernelWidth * SIMD_VECTOR_SIZE);
	const size_t colRemainder = a->cols % SIMD_VECTOR_SIZE;
	ThreadPool *pool = ThreadPool::getInstance();
	std::vector<TaskHandle> futures(nRows);
	// Do the main block
	//#pragma omp parallel num_threads(3)
	// OpenMP does not really seem to get much of a speedup out of the box. Is it too conservative with shared resources or something?
	for (int row = 0; row < nRows; row++) {
		auto rowTask = [=, &result]() {
			for (size_t col = 0; col < nCols; col++) {
				float *resultLocation = MatrixGetAddr(&result, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE);
				Microkernel(a, b, resultLocation, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, kernelHeight);
				if (colRemainder) {
					MicrokernelRemainder(a, b, resultLocation, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, kernelHeight);
				}
			}
			return true;
		};
		futures[row] = pool->spawnTask(rowTask);
	}

	for (size_t task = 0; task < nRows; task++) {
		pool->activeWait(futures[task]);
	}

	/*
	for (int row = 0; row < nRows; row++)
	{
		for (size_t col = 0; col < nCols; col++)
		{
			float *resultLocation = MatrixGetAddr(&result, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE);
			Microkernel(a, b, resultLocation, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, kernelHeight);
			if (colRemainder) {
				MicrokernelRemainder(a, b, &result, row * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, kernelHeight);
			}
		}
	}*/
	// Do the remaining lines if any (we may still have some columns to fill in afterwards)
	if (rowRemainder) {
		for (size_t col = 0; col < nCols; col++) {
			float *resultLocation = MatrixGetAddr(&result, nRows * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE);
			Microkernel(a, b, resultLocation, nRows * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, rowRemainder);
			if (colRemainder) {
				MicrokernelRemainder(a, b, resultLocation, nRows * kernelHeight, col * kernelWidth * SIMD_VECTOR_SIZE, kernelWidth, rowRemainder);
			}
		}
	}
	// Clean up what is left (columns)
	if (nCols * kernelWidth * SIMD_VECTOR_SIZE < b->cols)
	{
		for (size_t i = 0; i < a->rows; i++)
		{
			const float *aRow = a->contents + i * a->cols;
			float *resultRow = result.contents + i * result.cols;
			for (size_t k = 0; k < b->rows; k++)
			{
				const float *bRow = b->contents + k * b->cols;
				const float Aik = aRow[k];
				for (size_t j = nCols * kernelWidth * SIMD_VECTOR_SIZE; j < b->cols; j++)
				{
					//todo(AION): same as above, this needs to be compiler dependant to be in line
					//resultRow[j] += Aik * bRow[j];
					resultRow[j] = fma(bRow[j], Aik, resultRow[j]);

				}
			}
		}
	}
	return result;
}


#ifdef __APPLE__
typedef struct benchmark_endpoint {
	struct timespec time;
} benchmark_endpoint;

#else 
typedef struct benchmark_endpoint {
	clock_t benchmark_clock;
} benchmark_endpoint;
#endif

typedef struct benchmark_results
{
	double gflops;
	double total_time;
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
	return TimeDiffInSeconds(start.time, end.time);
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
	printf("Total time: %lf\n", benchmark.total_time);
	printf("===========\n");
}


double GFlops(double nOps, double timeInSeconds)
{
	return nOps * 2e-9 / timeInSeconds;
}

benchmark_results TestRandomNaive(size_t size_i, size_t size_j, size_t size_k, int nTrials, const size_t kernelWidth, const size_t kernelHeight)
{
	const matrix_f32 a = RandomMatrixf32(size_i, size_k);
	const matrix_f32 b = RandomMatrixf32(size_k, size_j);
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
			PrintMatrixf32(&targetResult);
			printf("\n\n\n");
			Matrixf32MultiplyInplace(&testResult, -1);
			matrix_f32 diff = Matrixf32Sum(&targetResult, &testResult);
			PrintMatrixf32(&diff);
		}
		free(testResult.contents);
	}
	averageResult.gflops = GFlops(nOps * nTrials, averageResult.total_time);
	return averageResult;
}

benchmark_results TestRandom(size_t size_i, size_t size_j, size_t size_k, int nTrials, const size_t kernelWidth, const size_t kernelHeight)
{
	const matrix_f32 a = RandomMatrixf32(size_i, size_k);
	const matrix_f32 b = RandomMatrixf32(size_k, size_j);
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
		matrix_f32 testResult = Matrixf32MicrokernelMultiply(&a, &b, kernelWidth, kernelHeight);
		benchmark_endpoint end = GetBenchmarkEndpoint();
		double timeInSeconds = GetTimeDiffFromEndpoints(start, end);
		averageResult.total_time += timeInSeconds;
		averageResult.max_error = GetMaxError(&testResult, &targetResult);
		averageResult.result = testResult;
		averageResult.gflops += GFlops(nOps, timeInSeconds);
		if (nTrials == 1) {
			PrintMatrixf32(&targetResult);
			printf("\n\n\n");
			Matrixf32MultiplyInplace(&testResult, -1);
			matrix_f32 diff = Matrixf32Sum(&targetResult, &testResult);
			PrintMatrixf32(&diff);
		}
		free(testResult.contents);
	}
	averageResult.gflops = GFlops(nOps * nTrials, averageResult.total_time);
	return averageResult;
}

benchmark_results TestOnes(size_t size_i, size_t size_j, size_t size_k, int nTrials, const size_t kernelWidth, const size_t kernelHeight)
{
	benchmark_results averageResult;
	const matrix_f32 a = Onesf32(size_i, size_k);
	const matrix_f32 b = Onesf32(size_k, size_j);
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
		matrix_f32 testResult = Matrixf32MicrokernelMultiply(&a, &b, kernelWidth, kernelHeight);
		benchmark_endpoint end = GetBenchmarkEndpoint();
		double timeInSeconds = GetTimeDiffFromEndpoints(start, end);
		averageResult.total_time += timeInSeconds;
		averageResult.max_error = GetMaxError(&testResult, &targetResult);
		averageResult.result = testResult;
		if (nTrials == 1) {
			PrintMatrixf32(&targetResult);
			printf("\n\n\n");
			Matrixf32MultiplyInplace(&testResult, -1);
			matrix_f32 diff = Matrixf32Sum(&targetResult, &testResult);
			PrintMatrixf32(&diff);
		}
		free(testResult.contents);
	}
	averageResult.gflops = GFlops(nOps * nTrials, averageResult.total_time);
	return averageResult;
}

int MatrixBenchmarksMain(size_t sizeI, size_t sizeJ, size_t sizeK, int nTrials, const size_t kernelWidth, const size_t kernelHeight)
{
	printf("Testing Naive multiply on random matrices\n");
	benchmark_results benchmark = TestRandomNaive(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight);
	PrintBenchmark(benchmark);
	printf("Testing microkernel multiply on random matrices\n");
	benchmark = TestRandom(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight);
	PrintBenchmark(benchmark);
	printf("Testing microkernel multiply on matrices of ones\n");
	benchmark = TestOnes(sizeI, sizeJ, sizeK, nTrials, kernelWidth, kernelHeight);
	PrintBenchmark(benchmark);
	return 0;
}
