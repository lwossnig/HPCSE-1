// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>
#include <iomanip>
#include <omp.h>

int main()
{
  unsigned long const nterms = 100000000;
  long double sum=0.;

  #pragma omp parallel reduction(+:sum)
  {
    int i = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    
    long double const step = (nterms+0.5l) /  nthreads;
    int j = (i+1) * step;
    for (std::size_t t = i * step; t < j; ++t)
      sum += (1.0 - 2* (t % 2)) / (2*t + 1);
  }
  
  std::cout << "pi=" << std::setprecision(18) << 4.*sum << std::endl;
  return 0;
}