#include <chrono>
#include "timer.hpp"

Timer::Timer(void)
	    {
		    begin = std::chrono::high_resolution_clock::now();
		    end = std::chrono::high_resolution_clock::now();
	    }
Timer::~Timer(void)
{}

void Timer::start()
	    {
		    begin = std::chrono::high_resolution_clock::now();
	    }
void Timer::stop()
	    {
		    end = std::chrono::high_resolution_clock::now();
	    }

double Timer::duration() // returns seconds of duration
	    {
		    std::chrono::duration< double > time = end - begin;
		    return time.count(); 
		    // returns time in seconds (time class-chrono duration!)
	    }

