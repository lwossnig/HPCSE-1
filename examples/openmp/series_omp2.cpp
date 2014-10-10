// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iomanip>
#include <iostream>

int main()
{
  unsigned long const nterms = 100000000;
  
  long double sum=0.;
  
  #pragma omp parallel shared(sum)
  {
    #pragma omp for reduction(+:sum)
    for (std::size_t t = 0; t < nterms; ++t)
      sum += (1.0 - 2* (t % 2)) / (2*t + 1);
  }
  
  std::cout << "pi=" << std::setprecision(18) << 4.*sum << std::endl;
  return 0;
}