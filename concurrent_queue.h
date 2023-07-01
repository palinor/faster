#pragma once
#include <queue>
#include <mutex>
#include <vector>

template <class T>
class ConcurrentQueue {
	std::queue<T> my_queue_;
	mutable std::mutex my_mutex_;
	std::condition_variable my_cv_;
	bool my_interrupt_;

public:

	ConcurrentQueue() : my_interrupt_(false) {}
	~ConcurrentQueue() { interrupt(); }

	void interrupt() {
		{
			std::lock_guard<std::mutex> lk(my_mutex_);
			my_interrupt_ = true;
		}
		my_cv_.notify_all();
	}

	void resetInterrupt() {
		my_interrupt_ = false;
	}

	bool empty() const {
		std::lock_guard<std::mutex> lk(my_mutex_);
		return my_queue_.empty();
	}

	void push(T t) {
		std::lock_guard<std::mutex> lk(my_mutex_);
		my_queue_.push(std::move(t));
		my_cv_.notify_one();
	}

	bool pop(T &t) {
		std::unique_lock<std::mutex> lk(my_mutex_);
		while (!my_interrupt_ && my_queue_.empty()) my_cv_.wait(lk);
		if (my_interrupt_) return false;
		t = std::move(my_queue_.front());
		my_queue_.pop();
		return true;
	}

	bool tryPop(T &t) {
		std::lock_guard<std::mutex> lk(my_mutex_);
		if (my_queue_.empty()) return false;
		t = std::move(my_queue_.front());
		my_queue_.pop();
		return true;
	}

	void clear() {
		std::queue<T> empty;
		swap(my_queue_, empty);
	}
};
