// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include "simpson.hpp"
#include <cmath>
#include <iostream>
#include <thread>
#include <future>

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
  
  // even easier: launch an asynchronous function call
  std::future<double> fi = std::async(simpson,func,a,a+(b-a)/2.,nsteps/2);
  
  // locally integrate the second half
  double result = simpson(func,a+(b-a)/2.,b,nsteps/2);
  
  std::cout << result+fi.get() << std::endl;
  
  return 0;
}