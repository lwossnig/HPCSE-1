// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include "simpson.hpp"
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
  
  // print the result
  std::cout << simpson(func,a,b,nsteps) << std::endl;
  
  return 0;
}