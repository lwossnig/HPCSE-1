// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>
#include <chrono>

int fibonacci(int n)
{
  if (n < 3) return 1;
  return fibonacci(n-1) + fibonacci(n-2);
}

int main()
{
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  start = std::chrono::high_resolution_clock::now();
  int result = fibonacci(42);
  end = std::chrono::high_resolution_clock::now();
  
  int elapsed_microseconds = std::chrono::duration_cast<std::chrono::microseconds>
  (end-start).count();
  
  std::cout << "elapsed time: " << elapsed_microseconds << "ms\n";
}

