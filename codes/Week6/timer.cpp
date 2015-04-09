#include "timer.hpp"
#include <chrono>

Timer::Timer() {
  tstart = std::chrono::high_resolution_clock::now();
}

void Timer::start() {
  tstart = std::chrono::high_resolution_clock::now();
}

void Timer::stop() {
  tend = std::chrono::high_resolution_clock::now();
}

double Timer::duration() const {
    
  return std::chrono::duration< double >(tend - tstart).count();
}
