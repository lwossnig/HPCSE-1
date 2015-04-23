#include <iostream>
#include <chrono>
#ifndef _TIMER_CLASS
#define _TIMER_CLASS

class Timer
{
	private:
		std::chrono::time_point< std::chrono::high_resolution_clock > begin , end;
	
	public:
		Timer(void);
		// constructor
		
		~Timer(); 
		//destructor

		void start();
		void stop();
		double duration(); // returns seconds of duration

};
#endif
