#include "simpson.hpp"
#include <iostream>
#include <cmath>

// a function with two variables
double expax(double a, double x)
{
  return std::exp(a*x);
}

int main()
{
  // where do we set a?
  std::cout << simpson(expax,0.,1.,100) << std::endl;
  return 0;
}

