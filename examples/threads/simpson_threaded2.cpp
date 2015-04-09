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
  
  // create a packaged task
  std::packaged_task<double()> pt(std::bind(simpson,func,a,a+(b-a)/2.,nsteps/2));
  std::future<double> fi = pt.get_future(); // get the future return value
  std::thread t (std::move(pt));            // launch the thread
  
  // locally integrate the second half
  double result = simpson(func,a+(b-a)/2.,b,nsteps/2);

  // wait for the task to finish and the future to be ready
  // fi.wait(); // not needed since it is implicit in fi.get() below
  
  std::cout << result+fi.get() << std::endl;
  t.join();
  
  return 0;
}