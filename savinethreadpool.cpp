#include "savinethreadpool.h"

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
		}
		else {
			f.wait();
		}
	}
	return i_did_work;
}


ThreadPool ThreadPool::my_instance_;
thread_local size_t ThreadPool::my_thread_number_ = 0;
