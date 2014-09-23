// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>
#include <thread>
#include <vector>

void printer( int n )
{
  for ( int i = 0; i < 100; ++i)
    std::cout << "do not garble thread " << n << ": " << i << std::endl;
}


int main()
{
  std::vector<std::thread> threads;
  
  for (int n = 1; n < 10; ++n)
    threads.push_back(std::thread(printer, n));
  
  for (std::thread& t : threads)
    t.join();
}

