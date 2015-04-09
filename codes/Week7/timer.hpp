#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
    
class Timer {
public:
  Timer();
  void start();
  void stop();
  double duration() const;

private:
  typedef std::chrono::high_resolution_clock::time_point time_point_t;
  time_point_t tstart, tend;
};

#endif // !defined TIMER_HPP
