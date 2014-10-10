// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>
#include <chrono>
#include <ctime>

int fibonacci(int n)
{
  if (n < 3) return 1;
  return fibonacci(n-1) + fibonacci(n-2);
}

int main()
{
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  int result = fibonacci(42);
  end = std::chrono::system_clock::now();
  
  int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
  (end-start).count();
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  
  std::cout << "finished computation at " << std::ctime(&end_time)
  << "elapsed time: " << elapsed_seconds << "ms\n";
}

