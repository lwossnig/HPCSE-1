#include <iostream>
#include <omp.h>

int main()
{
  #pragma omp parallel
  {
    // now we execute this block in multiple threads
    #pragma omp critical
    std::cout << "I am thread " << omp_get_thread_num()
              << " of " << omp_get_num_threads() << " threads." << std::endl;
  }
  
}