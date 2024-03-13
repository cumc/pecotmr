#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <cassert>
#include <thread>

class Function_pool {
    private:
	std::vector<std::thread> thread_pool;
	std::queue<std::function<void()>> m_function_queue;
	std::mutex m_lock;
	std::condition_variable m_data_condition;
	std::condition_variable finished;
	std::atomic<bool> m_accept_functions;
	int busy;
	bool stop;

    public:
	Function_pool(unsigned int n);
	~Function_pool();
	void push(std::function<void()> func);
	void infinite_loop_func();
	void waitFinished();
};
