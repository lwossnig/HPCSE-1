#include <random>
#include <iostream>


int main()
{
  // create an engine
  std::mt19937 mt;
  
  // seed the generator
  mt.seed(42);
  
  // create four distributions
  std::uniform_int_distribution<int>     uint_d(0,10);
  std::uniform_real_distribution<double> ureal_d(0.,10.);
  std::normal_distribution<double>       normal_d(0.,4.);
  std::exponential_distribution<double>  exp_d(1.);
  
  // create random numbers with each of these distributions:
  
  std::cout << uint_d(mt)   << "\n";
  std::cout << ureal_d(mt)  << "\n";
  std::cout << normal_d(mt) << "\n";
  std::cout << exp_d(mt)    << "\n";
    
  return 0;
}