#include "simpson.hpp"
#include <iostream>
#include <cmath>

// an ugly global variable
double a;

// the function to be integrated
double expax(double x)
{
  return std::exp(a*x);
}

int main()
{
  a=3.4;
  std::cout << simpson(expax,0.,1.,100) << std::endl;
  return 0;
}

