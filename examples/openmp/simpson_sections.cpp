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

  #pragma omp parallel shared(result)
  {
    #pragma omp sections reduction(+:result)
    {
      #pragma omp section 
      {
        result = simpson(func,a,a+0.5*(b-a),nsteps/2);
      }
      #pragma omp section
      {
        result = simpson(func,a+0.5*(b-a),b,nsteps/2);
      }
    }
  }
  
  std::cout << result << std::endl;
  
  return 0;
}