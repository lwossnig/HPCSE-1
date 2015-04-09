// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <vector>
#include <iostream>
#include <thread>
#include <iomanip>


// sum terms [i-j) of the power series for pi/4
void sumterms(long double& sum, std::size_t i, std::size_t j)
{  
  for (std::size_t t = i; t < j; ++t)
    sum += (1.0 - 2* (t % 2)) / (2*t + 1);
}


int main()
{
  // decide how many threads to use
  std::size_t const nthreads = std::max(1u, std::thread::hardware_concurrency());

  std::vector<std::thread> threads(nthreads);
  // let us just use a single result
  long double result=0.;

  unsigned long const nterms = 100000000;
  long double const step = (nterms+0.5l) /  nthreads;
  
  for (unsigned i = 0; i < nthreads; ++i)
    threads[i] = std::thread(sumterms,std::ref(result), i * step, (i+1) * step);
  
  for (std::thread& t : threads)
    t.join();
    
  std::cout << "pi=" << std::setprecision(18) << 4.*result << std::endl;
  return 0;
}