#include <chrono>
#include <condition_variable>
#include <future>
#include <queue>
#include <mutex>
#include <thread>
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
