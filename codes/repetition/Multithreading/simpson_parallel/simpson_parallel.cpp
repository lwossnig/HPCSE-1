/// simpson.cpp
#include "simpson_parallel.hpp"
#include <cassert>
#include "sync_io.hpp"

void simpson(const double a, const double b,const unsigned bins, std::function<double(double)> f, std::pair<double, std::mutex>& result)
{
	assert(bins > 0);
	assert(f!=NULL);

	sync(std::cout) << "Thread starts with " << a << " and stops with " << b << std::endl;

	const unsigned int steps = 2*bins+1;

	const double dr = (b-1) / (steps-1);
	double I = f(a);

	for(unsigned int i=1; i < steps-1; ++i){
		I += 2* (1.0 + i%2) * f(a+dr*i);
	}

	I += f(b);
	I *= (1./3.) * dr;

	std::lock_guard<std::mutex> l(result.second);
	result.first += I;
}

