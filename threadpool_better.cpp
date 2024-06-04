#include "threadpool_better.h"
// This should all be doable without going through futures and promises
// How much of the cpp features can I remove from my multithreading and still have it work?

void concurrentTaskQueueInit(ConcurrentTaskQueue *queue, size_t max_size) {
	queue->tasks = (TaskInfoWrapper *)malloc(max_size * sizeof(TaskInfoWrapper));
	queue->max_size = max_size;
	queue->starting_max_size = max_size;
	queue->current_size = 0;
	queue->current_head = 0;
	queue->is_interrupted = false;
}

void concurrentTaskQueueReset(ConcurrentTaskQueue *queue) {
	free(queue->tasks);
	queue->tasks = (TaskInfoWrapper *)malloc((queue->starting_max_size) * sizeof(TaskInfoWrapper));
	queue->max_size = queue->starting_max_size;
	queue->current_size = 0;
	queue->current_head = 0;
	queue->is_interrupted = false;
}


void concurrentTaskQueuePush(void *task_info, ConcurrentTaskQueueCallback *task_function, ConcurrentTaskQueue *queue) {
	queue->mutex.lock();
	if (queue->current_size == queue->max_size) {
		size_t new_size = queue->max_size * 2 + 1;
		TaskInfoWrapper *new_tasks = reinterpret_cast<TaskInfoWrapper *>(realloc(queue->tasks, new_size * sizeof(TaskInfoWrapper)));
		if (!new_tasks) {
			exit(1);
		}
		queue->tasks = new_tasks;
		queue->max_size = new_size;
	}
	TaskInfoWrapper *next_task = queue->tasks + queue->current_size;
	next_task->task_info = task_info;
	next_task->task_function = task_function;
	queue->current_size++;
	queue->cv.notify_one();
	queue->mutex.unlock();
}

inline bool concurrentTaskQueueIsEmpty(ConcurrentTaskQueue *queue) {
	return (queue->current_head == queue->current_size);
}

volatile TaskInfoWrapper *concurrentTaskQueuePopTask(ConcurrentTaskQueue *queue) {
	std::unique_lock<std::mutex> lk(queue->mutex);
	while (!(queue->is_interrupted) && concurrentTaskQueueIsEmpty(queue)) {
		queue->cv.wait(lk);
	}
	if ((queue->is_interrupted)) {
		return nullptr;
	}
	if (concurrentTaskQueueIsEmpty(queue)) {
		return nullptr;
	}
	volatile TaskInfoWrapper *return_task = queue->tasks + queue->current_head;
	queue->current_head++;
	return return_task;
}

void concurrentTaskQueueInterrupt(ConcurrentTaskQueue *queue) {
	queue->mutex.lock();
	queue->is_interrupted = true;
	queue->mutex.unlock();
	queue->cv.notify_all();
}


TaskInfoWrapper *concurrentTaskQueueTryPop(ConcurrentTaskQueue *queue) {
	queue->mutex.lock();
	if (concurrentTaskQueueIsEmpty(queue)) {
		queue->mutex.unlock();
		return nullptr;
	}
	TaskInfoWrapper *result_task = queue->tasks + queue->current_head++;
	queue->mutex.unlock();
	return result_task;
}



void threadPoolRunThreadFunc(ThreadPool *pool) {
	while (!(pool->is_interrupted)) {
		volatile TaskInfoWrapper *task_wrapper = concurrentTaskQueueTryPop(&(pool->task_queue));
		if (!(pool->is_interrupted) && !!(task_wrapper)) {
			task_wrapper->task_function(task_wrapper->task_info);
		}
	}
}

void threadPoolInit(ThreadPool *pool, size_t starting_queue_size, size_t n_threads) {
	pool->is_interrupted = false;
	if (pool->threads) {
		free(pool->threads);
	}
	pool->threads = (std::thread *)malloc(n_threads * sizeof(std::thread));
	pool->number_of_total_threads = n_threads;
	for (size_t thread_number = 0; thread_number < n_threads; thread_number++) {
		new(pool->threads + thread_number) std::thread(threadPoolRunThreadFunc, pool);
	}
	concurrentTaskQueueInit(&(pool->task_queue), starting_queue_size);
	pool->is_active = true;
}

void threadPoolInterrupt(ThreadPool *pool) {
	if (pool->is_active) {
		pool->is_interrupted = true;
		concurrentTaskQueueInterrupt(&(pool->task_queue));
	}
	for (size_t thread_number = 0; thread_number < pool->number_of_total_threads; thread_number++) {
		pool->threads[thread_number].join();
	}
	free(pool->threads);
	pool->is_active = false;
	pool->is_interrupted = false;
	pool->task_queue.is_interrupted = false;
	if (pool->task_queue.tasks) {
		free(pool->task_queue.tasks);
	}
}

bool threadPoolActiveWait(ThreadPool *pool) {
	bool i_did_work = false;
	while (!concurrentTaskQueueIsEmpty(&(pool->task_queue))) {
		if (volatile TaskInfoWrapper *task_wrapper = concurrentTaskQueueTryPop(&(pool->task_queue))) {
			task_wrapper->task_function(task_wrapper->task_info);
			i_did_work = true;
		}
	}
	return i_did_work;
}

