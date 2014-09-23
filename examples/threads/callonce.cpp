// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>
#include <thread>
#include <vector>

std::once_flag printonce_flag;

void printonce()
{
  std::cout << "This should be printed only once\n";
}


int main()
{
  std::vector<std::thread> threads;
  
  for (int n = 0; n < 10; ++n)
    threads.push_back(
        std::thread([&](){std::call_once(printonce_flag,printonce);}));
  
  for (std::thread& t : threads)
    t.join();

  // and again
  for (int n = 0; n < 10; ++n)
    threads[n] =
        std::thread([&](){std::call_once(printonce_flag,printonce);});

  for (std::thread& t : threads)
    t.join();
  
  // and again
  std::call_once(printonce_flag,printonce);
}

