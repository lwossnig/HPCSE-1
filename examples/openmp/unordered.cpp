// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>

int main()
{
  #pragma omp parallel for
  for (int i=0; i < 100; ++i) {
    // do some (fake) work
    int j=i;
    #pragma omp critical
    std::cout << "Hello from the " << j << "-th iteration\n";
  }
}

    
 