#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>
#include "timer.hpp"

double f(double x)
{
	return std::log(x)*std::sqrt(x);
}

// WolframAlpha: integral_1^4 log(x) sqrt(x) dx = 4/9 (4 log(64)-7) ~~ 4.28245881486164
int main(int argc, char *argv[])
{
	double a = 1.0;
	double b = 4.0;
	unsigned long const n = 200UL*35389440UL;
	const double dx = (b-a)/n;

	unsigned int nthreads = 1;
	if (argc > 1) nthreads = atoi(argv[1]);
	unsigned long nsteps = n / nthreads;

	std::vector<std::thread> threads(nthreads);

	double S = 0;
	std::mutex S_mutex;			// mutual exclusion

	timer t;
	t.start();
	for (unsigned thr = 0; thr < nthreads; thr++) {
#if DBG
		std::cout << "spawning thread " << thr << std::endl;
#endif
		threads[thr] = std::thread([&,thr]() {
			double local_S = 0;
			double x0 = a + (thr * nsteps + 0.5)*dx;
			for (unsigned long i = 0; i < nsteps; i++) {
				double xi = x0 + i*dx; 
				local_S += f(xi);
			}
			local_S *= dx;
			S_mutex.lock();
			S += local_S;
			S_mutex.unlock();
		});
	}

	for (std::thread& thr : threads)
		thr.join();

	t.stop();

	std::cout << "Threads=" << nthreads << " Time=" << t.get_timing() << " seconds, Result=" << std::setprecision(8) << S << std::endl;

	return 0; 
}
