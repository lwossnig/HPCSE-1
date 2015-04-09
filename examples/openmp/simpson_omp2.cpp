// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include "simpson.hpp"
#include <omp.h>
#include <cmath>
#include <iostream>

// The function to integrate
double func(double x)
{
  return x * std::sin(x);
}

int main()
{
  double a;            // lower bound of integration
  double b;            // upper bound of integration
  unsigned int nsteps; // number of subintervals for integration
  
  // read the parameters
  std::cin >> a >> b >> nsteps;
  
  double result=0.;

  #pragma omp parallel
  {
    int i = omp_get_thread_num();
    int n = omp_get_num_threads();
    double delta = (b-a)/n;
    // integrate just one part in each thread
    double r = simpson(func,a+i*delta,a+(i+1)*delta,nsteps/n);
    #pragma omp critical (simpsonresult)
    result += r;
  }
  
  std::cout << result << std::endl;
  
  return 0;
}