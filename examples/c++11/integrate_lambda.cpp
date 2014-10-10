#include "simpson.hpp"
#include <iostream>
#include <cmath>

int main()
{
  double a=3.4;
  
  // create a lambda function
  // [=] indicates that the variable a should be used inside the lambda
  auto f = [=] (double x) { return std::exp(a*x); };
  
  std::cout << simpson(f,0.,1.,100) << std::endl;
  return 0;
}

