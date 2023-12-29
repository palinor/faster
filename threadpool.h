#pragma once
#include <mutex>
#include <thread>

using namespace std::chrono_literals;

typedef void concurrent_taskqueue_callback(void *taskInfo);

struct taskinfo_wrapper {
	void *taskInfo_;
	concurrent_taskqueue_callback *taskFunction_;
};

struct concurrent_taskinfo_queue {
	taskinfo_wrapper *tasks_;
	size_t volatile startingMaxSize_;
	size_t volatile maxSize_;
	size_t volatile currentSize_;
	size_t volatile currentHead_;
	std::mutex mutex_;
	std::condition_variable cv_;
	bool interrupt_;
};

struct thread_pool {
	std::thread *threads_ = nullptr;
	volatile size_t nTotalThreads_ = 0;
	volatile bool isActive_ = false;
	volatile bool isInterrupt_ = false;
	concurrent_taskinfo_queue queue_;
};

void InitThreadPool(
	thread_pool *input, 
	size_t initQueueSize,
	size_t nThreads = std::thread::hardware_concurrency() - 1
);
void RunThreadFunc(thread_pool *pool);
void StopThreadPool(thread_pool *pool);
bool ThreadPoolActiveWait(thread_pool *pool);
void InitConcurrentTaskInfoQueue(concurrent_taskinfo_queue *queue, size_t maxSize);
void ResetConcurrentTaskQueue(concurrent_taskinfo_queue *queue);
void PushTaskToQueue(void *taskInfo, concurrent_taskqueue_callback *taskFunction, concurrent_taskinfo_queue *queue);
inline bool IsQueueEmpty(concurrent_taskinfo_queue *queue);
volatile taskinfo_wrapper *PopTaskFromQueue(concurrent_taskinfo_queue *queue);
void InterruptConcurrentTaskInfoQueue(concurrent_taskinfo_queue *input);
taskinfo_wrapper *TryPopTaskFromQueue(concurrent_taskinfo_queue *queue);
void RunThreadFunc(thread_pool *pool);
void InitThreadPool(thread_pool *pool, size_t initQueueSize, size_t nThreads);
void StopThreadPool(thread_pool *pool);
bool ThreadPoolActiveWait(thread_pool *pool);

