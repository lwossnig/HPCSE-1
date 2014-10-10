#include "simpson.hpp"
#include <iostream>
#include <cmath>
#include <functional>

// a function with two variables
double expax(double a, double x)
{
  return std::exp(a*x);
}

int main()
{
  using namespace std::placeholders;
  
  double a=3.4;
  // bind one argument
  // _1, _2, .... are used for unbound arguments of the resulting function
  auto f = std::bind(expax,3.4,_1);
  std::cout << simpson(f,0.,1.,100) << std::endl;
  return 0;
}

