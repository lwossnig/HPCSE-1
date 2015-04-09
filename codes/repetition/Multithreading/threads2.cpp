#include <iostream>
#include <vector>
#include <numeric>
#include <random>
#include <thread>
#include <mutex>
#include "sync_io.hpp"

void sumterms(std::pair<double, std::mutex>& result, std::size_t i, std::size_t j)
{
	sync(std::cout) << "thread[" << i << "] working now!" << std::endl;
	double sum = 0.0;
	for(std::size_t t=i; t<j;++t){
		sum += (1.0 - 2* (t % 2)) / (2 * t + 1);
	}
	std::lock_guard<std::mutex> l(result.second);
	result.first+= sum;
}

int main()
{
	std::size_t const nthreads = std::max(1u, std::thread::hardware_concurrency());
	std::pair<double, std::mutex> sum;
	std::vector<std::thread> threads(nthreads);
	for (unsigned i=0; i < threads.size(); i++){
		threads[i] = std::thread(
				sumterms, std::ref(sum), 
				i*1000, (i+1)*1000);
	}
	for(auto& t : threads)
		t.join();
	
	std::cout << sum.first << std::endl;
}

