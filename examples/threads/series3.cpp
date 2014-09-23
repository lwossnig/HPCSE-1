// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <vector>
#include <iostream>
#include <thread>
#include <iomanip>
#include <mutex>

// sum terms [i-j) of the power series for pi/4
void sumterms(std::pair<long double, std::mutex>& result, std::size_t i, std::size_t j)
{
  long double sum=0.;
  
  for (std::size_t t = i; t < j; ++t)
    sum += (1.0 - 2* (t % 2)) / (2*t + 1);
  
  std::lock_guard<std::mutex> l (result.second);
  result.first += sum;
}


int main()
{
  // decide how many threads to use
  std::size_t const nthreads = std::max(1u, std::thread::hardware_concurrency());

  std::vector<std::thread> threads(nthreads);
  // let us just use a single result
  std::pair<long double, std::mutex> result;
  result.first = 0.;
  
  unsigned long const nterms = 100000000;
  long double const step = (nterms+0.5l) /  nthreads;
  
  for (unsigned i = 0; i < nthreads; ++i)
    threads[i] = std::thread(sumterms,std::ref(result), i * step, (i+1) * step);
  
  for (std::thread& t : threads)
    t.join();
    
  // no need to lock here
  std::cout << "pi=" << std::setprecision(18) << 4.*result.first << std::endl;
  return 0;
}