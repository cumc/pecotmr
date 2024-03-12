
#include "function_pool.h"
#include <iostream>
/* based on 
https://stackoverflow.com/a/51400041

https://stackoverflow.com/questions/23896421/efficiently-waiting-for-all-tasks-in-a-threadpool-to-finish
*/

Function_pool::Function_pool(unsigned int n) : m_function_queue(), m_lock(), m_data_condition(), busy(0), stop(false) {
    for (size_t i=0; i<n; i++) {
	thread_pool.push_back(std::thread(&Function_pool::infinite_loop_func, this));
    }
}

Function_pool::~Function_pool() {
    std::unique_lock<std::mutex> lock(m_lock);
    stop = true;
    m_data_condition.notify_all();
    lock.unlock();
    for (size_t i=0; i<thread_pool.size(); i++) {
	thread_pool[i].join();
    }
}

void Function_pool::push(std::function<void()> func) {
    std::unique_lock<std::mutex> lock(m_lock);
    m_function_queue.push(func);
    lock.unlock();
    m_data_condition.notify_one();
}

void Function_pool::infinite_loop_func() {
    std::function<void()> func;
    while (true) {
	std::unique_lock<std::mutex> lock(m_lock);
	m_data_condition.wait(lock, [this]() {return \
		!m_function_queue.empty() || stop; });
	if (!m_function_queue.empty()) {
	    busy++;
	    func = m_function_queue.front();
	    m_function_queue.pop();
	    lock.unlock();
	    func();
	    lock.lock();
	    busy--;
	    lock.unlock();
	    finished.notify_one();
	}
	else if (stop) {
	    return;
	}
    }
}

void Function_pool::waitFinished() {
    std::unique_lock<std::mutex> lock(m_lock);
    finished.wait(lock, [this](){return m_function_queue.empty() && busy == 0; });
}

