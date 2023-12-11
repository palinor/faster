#pragma once
#include <chrono>
#include <condition_variable>
#include <future>
#include <queue>
#include <mutex>
#include <thread>
#include <vector>

#include "savine_concurrent_queue.h"

using Task = std::packaged_task<bool(void)>;
using TaskHandle = std::future<bool>;

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

