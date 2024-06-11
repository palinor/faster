#pragma once
#include <condition_variable>
#include <mutex>
#include <thread>
#include "platform_build.h"

using namespace std::chrono_literals;

typedef void ConcurrentTaskQueueCallback(void *task_info);

struct TaskInfoWrapper {
	void *task_info;
	ConcurrentTaskQueueCallback *task_function;
};

struct ConcurrentTaskQueue {
	TaskInfoWrapper *tasks;
	bool *task_is_done;
	size_t volatile starting_max_size;
	size_t volatile max_size;
	size_t volatile current_size;
	size_t volatile current_head;
	std::mutex mutex;
	std::condition_variable cv;
	bool is_interrupted;
};

struct ThreadPool {
	std::thread *threads = nullptr;
	volatile size_t number_of_total_threads = 0;
	volatile bool is_active = false;
	volatile bool is_interrupted = false;
	ConcurrentTaskQueue task_queue;
};

void threadPoolInit(
	ThreadPool *input, 
	size_t starting_queue_size,
	size_t number_of_threads = std::thread::hardware_concurrency() - 1
);
void threadPoolRunThreadFunc(ThreadPool *pool);
bool threadPoolActiveWait(ThreadPool *pool);
void concurrentTaskQueueInit(ConcurrentTaskQueue *queue, size_t max_size);
void concurrentTaskQueueReset(ConcurrentTaskQueue *queue);
void concurrentTaskQueuePush(void *task_info, ConcurrentTaskQueueCallback *task_function, ConcurrentTaskQueue *queue);
inline bool concurrentTaskQueueIsEmpty(ConcurrentTaskQueue *queue);
volatile TaskInfoWrapper *concurrentTaskQueuePopTask(ConcurrentTaskQueue *queue);
void concurrentTaskQueueInterrupt(ConcurrentTaskQueue *input);
TaskInfoWrapper *concurrentTaskQueueTryPop(ConcurrentTaskQueue *queue);
void threadPoolRunThreadFunc(ThreadPool *pool);
void threadPoolInit(ThreadPool *pool, size_t initQueueSize, size_t nThreads);
void threadPoolInterrupt(ThreadPool *pool);
bool threadPoolActiveWait(ThreadPool *pool);

