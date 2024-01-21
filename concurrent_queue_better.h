#pragma once
#include <future>
#include <functional>
#include <mutex>

using task = std::packaged_task<bool(void)>;
using task_handle = std::future<bool>;

struct concurrent_task_queue {
	task *tasks_;
	size_t startingMaxSize_;
	size_t maxSize_;
	size_t currentSize_;
	size_t currentHead_;
	mutable std::mutex mutex_;
	std::condition_variable cv_;
	bool interrupt_;
};

void InitConcurrentTaskQueue(concurrent_task_queue *queue, size_t startingMaxSize);
void InterruptConcurrentTaskQueue(concurrent_task_queue *queue);
inline bool IsQueueEmpty(concurrent_task_queue *queue);
task_handle PushTaskToQueue(std::function<bool(void)> taskFunction, concurrent_task_queue *queue);
bool PopTaskFromQueue(task *resultTask, concurrent_task_queue *queue);
bool TryPopTaskFromQueue(task *resultTask, concurrent_task_queue *queue);
void ResetConcurrentTaskQueue(concurrent_task_queue *queue);
