// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>
#include <mutex>
#include "omp_mutex.hpp"

int main()
{
  
  omp_mutex m;
  
  #pragma omp parallel for
  for (int i=0; i < 100; ++i) {
    {
      std::lock_guard<omp_mutex> lock(m);
      std::cout << "Hello from the " << i << "-th iteration\n";
    }
  }
}

    
 