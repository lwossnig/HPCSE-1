#include <random>
#include <iostream>

template <class Engine, class Distribution>
void print_random_numnbers(Engine& eng, Distribution& dist, std::string doc)
{
  std::cout << "\n" << doc << "\n";
  for (int i=0; i<10; ++i)
    std::cout << dist(eng) << "\n";
}


int main()
{
  // create two engines
  std::minstd_rand lcg;
  std::mt19937 mt;
  
  // seed the Mersenne twister using the lcg.
  // This can be used for stochastic seeding of multiple generators in parallel
  mt.seed(lcg());
  
  std::cout << lcg() << " " << mt() << std::endl;
  
  // create various distributions
  std::uniform_int_distribution<int>     uint_d(0,10);
  std::uniform_real_distribution<double> ureal_d(0.,10.);
  std::normal_distribution<double>       normal_d(0.,4.);
  std::exponential_distribution<double>  exp_d(1.);
  
  print_random_numnbers(mt,uint_d,"uniform int");
  print_random_numnbers(mt,ureal_d,"uniform real");
  print_random_numnbers(mt,normal_d,"normal");
  print_random_numnbers(mt,exp_d,"exponential");
  
  return 0;
}