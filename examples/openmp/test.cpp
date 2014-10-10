// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>
#include "mutex.hpp"

int main()
{
  
  mutex m;
  
  #pragma omp parallel for
  for (int i=0; i < 100; ++i) {
    {
      lock_guard<mutex> lock(m);
      std::cout << "Hello from the " << i << "-th iteration\n";
    }
  }
}

    
 
