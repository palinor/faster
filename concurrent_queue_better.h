#pragma once
#include <future>
#include <functional>
#include <mutex>

using Task = std::packaged_task<bool(void)>;
using TaskHandle = std::future<bool>;

struct ConcurrentTaskQueue {
	Task *tasks_;
	size_t startingMaxSize_;
	size_t maxSize_;
	size_t currentSize_;
	size_t currentHead_;
	mutable std::mutex mutex_;
	std::condition_variable cv_;
	bool interrupt_;
};

void InitConcurrentTaskQueue(ConcurrentTaskQueue *queue, size_t startingMaxSize);
void InterruptConcurrentTaskQueue(ConcurrentTaskQueue *queue);
inline bool IsQueueEmpty(ConcurrentTaskQueue *queue);
TaskHandle PushTaskToQueue(std::function<bool(void)> taskFunction, ConcurrentTaskQueue *queue);
bool PopTaskFromQueue(Task *resultTask, ConcurrentTaskQueue *queue);
bool TryPopTaskFromQueue(Task *resultTask, ConcurrentTaskQueue *queue);
void ResetConcurrentTaskQueue(ConcurrentTaskQueue *queue);
