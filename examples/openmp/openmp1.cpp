#include <iostream>
#include <omp.h>

int main()
{
  std::cout << "I am thread " << omp_get_thread_num()
            << " of " << omp_get_num_threads() << " threads." << std::endl;
}