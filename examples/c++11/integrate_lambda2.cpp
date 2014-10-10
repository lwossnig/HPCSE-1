#include "simpson.hpp"
#include <iostream>
#include <cmath>

int main()
{
  double a=3.4;
  
  std::cout << simpson([=] (double x) { return std::exp(a*x); },0.,1.,100)
            << std::endl;
  return 0;
}

