#include "simpson.hpp"
#include <iostream>
#include <cmath>

// a function object for exp(a*x)
class expax
{
public:
  // set the parameter a in the constructor
  expax(double a) : a_(a) {}
  
  // the function call operator calculates the function
  double operator()(double x) { return std::exp(a_*x);}
  
private:
  double a_; // the fixed parameter a
};

int main()
{
  double a=3.4;
  std::cout << simpson(expax(a),0.,1.,100) << std::endl;
  return 0;
}

