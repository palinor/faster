#include "concurrent_queue.h"

// This queue is not exactly memory safe. It needs to be manually reset every time we run a batch of jobs
void InitConcurrentTaskQueue(concurrent_task_queue *queue, size_t maxSize) {
	queue->tasks_ = (task *)malloc(maxSize * sizeof(task));
	queue->maxSize_ = maxSize;
	queue->startingMaxSize_ = maxSize;
	queue->currentSize_ = 0;
	queue->currentHead_ = 0;
	queue->interrupt_ = false;
}

void ResetConcurrentTaskQueue(concurrent_task_queue *queue) {
	free(queue->tasks_);
	queue->tasks_ = (task *)malloc((queue->startingMaxSize_) * sizeof(task));
	queue->maxSize_ = queue->startingMaxSize_;
	queue->currentSize_ = 0;
	queue->currentHead_ = 0;
	queue->interrupt_ = false;
}


task_handle PushTaskToQueue(std::function<bool(void)> taskFunction, concurrent_task_queue *queue) {
	queue->mutex_.lock();
	if (queue->currentSize_ == queue->maxSize_) {
		size_t newSize = queue->maxSize_ * 2 + 1;
		task *newTasks = reinterpret_cast<task *>(realloc(queue->tasks_, newSize * sizeof(task)));
		if (!newTasks) {
			exit(1);
		}
		queue->tasks_ = newTasks;
		queue->maxSize_ = newSize;
	}
	new(queue->tasks_ + (queue->currentSize_++)) task(std::move(taskFunction));
	task_handle handle = queue->tasks_[queue->currentSize_ - 1].get_future();
	//queue->tasks_[queue->currentSize_++] = task(std::move(element));
	queue->cv_.notify_one();
	queue->mutex_.unlock();
	return handle;
}

inline bool IsQueueEmpty(concurrent_task_queue *queue) {
	return (queue->currentHead_ == queue->currentSize_);
}

bool PopTaskFromQueue(task *outputTask, concurrent_task_queue *queue) {
	std::unique_lock<std::mutex> lk(queue->mutex_);
	while (!(queue->interrupt_) && IsQueueEmpty(queue)) {
		queue->cv_.wait(lk);
	}
	if ((queue->interrupt_)) {
		return false;
	}
	if (IsQueueEmpty(queue)) {
		return false;
	}
	*outputTask = std::move(queue->tasks_[queue->currentHead_++]);
	return true;
}

void InterruptConcurrentTaskQueue(concurrent_task_queue *input) {
	input->mutex_.lock();
	input->interrupt_ = true;
	input->mutex_.unlock();
	input->cv_.notify_all();
}


bool TryPopTaskFromQueue(task *resultTask, concurrent_task_queue *queue) {
	queue->mutex_.lock();
	if (IsQueueEmpty(queue)) {
		queue->mutex_.unlock();
		return false;
	}
	*resultTask = std::move(queue->tasks_[queue->currentHead_++]);
	queue->mutex_.unlock();
	return true;
}


