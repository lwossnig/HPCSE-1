// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>
#include <omp.h>

int main() {
  
  #pragma omp parallel
  #pragma omp single
  std::cout << "Only thread " << omp_get_thread_num()
            << " of " << omp_get_num_threads() << " is printing.\n";
  
  return 0;
}

