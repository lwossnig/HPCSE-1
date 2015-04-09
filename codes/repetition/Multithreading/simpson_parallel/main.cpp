#include <iostream>
#include <vector>
#include <thread>
#include <cmath>
#include <mutex>
#include "simpson_parallel.hpp"


double exp_minus_x(double x) {return std::exp(-x);}

int main()
{
	const unsigned nthreads = std::max(1u, std::thread::hardware_concurrency());
	std::pair<double, std::mutex> result;
	std::vector<std::thread> threads(nthreads);
	
	const unsigned bins = 100;
	const double a=0,b=10;
	double step_size =(double)(b-a)/nthreads;

	for(std::size_t i=0; i <nthreads; ++i){
		threads[i] = std::thread(simpson, i*step_size, (i+1)*step_size, bins, exp_minus_x, std::ref(result));
	}

	for(auto& t : threads){
		t.join();
	}

	std::cout << "Result of the integration from " << a << " till " << b << " of the chosen function (standard exp(-x) is: " << result.first << std::endl;

	return 0;
}

