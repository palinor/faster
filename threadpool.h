#pragma once
#include <chrono>
#include <functional>
#include <future>
#include <mutex>
#include <thread>


using namespace std::chrono_literals;


void InitThreadPool(
	thread_pool *input, 
	size_t initQueueSize,
	size_t nThreads = std::thread::hardware_concurrency() - 1
);
void RunThreadFunc(thread_pool *pool);
void StopThreadPool(thread_pool *pool);
bool ThreadPoolActiveWait(thread_pool *pool, const task_handle &f);
// inline task_handle PushTaskToThreadPool(std::function<bool(void)> taskFunction, thread_pool *pool);