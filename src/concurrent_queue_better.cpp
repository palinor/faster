#include "concurrent_queue_better.h"

// This queue is not exactly memory safe. It needs to be manually reset every time we run a batch of jobs
void InitConcurrentTaskQueue(ConcurrentTaskQueue *queue, size_t maxSize) {
	queue->tasks_ = (Task *)malloc(maxSize * sizeof(Task));
	queue->maxSize_ = maxSize;
	queue->startingMaxSize_ = maxSize;
	queue->currentSize_ = 0;
	queue->currentHead_ = 0;
	queue->interrupt_ = false;
}

void ResetConcurrentTaskQueue(ConcurrentTaskQueue *queue) {
	free(queue->tasks_);
	queue->tasks_ = (Task *)malloc((queue->startingMaxSize_) * sizeof(Task));
	queue->maxSize_ = queue->startingMaxSize_;
	queue->currentSize_ = 0;
	queue->currentHead_ = 0;
	queue->interrupt_ = false;
}


TaskHandle PushTaskToQueue(std::function<bool(void)> taskFunction, ConcurrentTaskQueue *queue) {
	queue->mutex_.lock();
	if (queue->currentSize_ == queue->maxSize_) {
		size_t newSize = queue->maxSize_ * 2 + 1;
		Task *newTasks = reinterpret_cast<Task *>(realloc(queue->tasks_, newSize * sizeof(Task)));
		if (!newTasks) {
			exit(1);
		}
		queue->tasks_ = newTasks;
		queue->maxSize_ = newSize;
	}
	new(queue->tasks_ + (queue->currentSize_++)) Task(std::move(taskFunction));
	TaskHandle handle = queue->tasks_[queue->currentSize_ - 1].get_future();
	//queue->tasks_[queue->currentSize_++] = task(std::move(element));
	queue->cv_.notify_one();
	queue->mutex_.unlock();
	return handle;
}

inline bool IsQueueEmpty(ConcurrentTaskQueue *queue) {
	return (queue->currentHead_ == queue->currentSize_);
}

bool PopTaskFromQueue(Task *outputTask, ConcurrentTaskQueue *queue) {
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

void InterruptConcurrentTaskQueue(ConcurrentTaskQueue *input) {
	input->mutex_.lock();
	input->interrupt_ = true;
	input->mutex_.unlock();
	input->cv_.notify_all();
}


bool TryPopTaskFromQueue(Task *resultTask, ConcurrentTaskQueue *queue) {
	queue->mutex_.lock();
	if (IsQueueEmpty(queue)) {
		queue->mutex_.unlock();
		return false;
	}
	*resultTask = std::move(queue->tasks_[queue->currentHead_++]);
	queue->mutex_.unlock();
	return true;
}


