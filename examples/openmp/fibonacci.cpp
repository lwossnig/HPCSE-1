// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>

int fibonacci(int n)
{
  int i, j;
  
  if (n<2)
    return n;
  else {
    #pragma omp task shared(i) firstprivate(n) untied final (n<=5)
    i = fibonacci(n-1);
    
    #pragma omp task shared(j) firstprivate(n) untied final (n<=5)
    j = fibonacci(n-2);
    
    #pragma omp taskwait
    return i + j;
  }
}

int main()
{
  int n;
  std::cin >> n;

  #pragma omp parallel shared(n)
  {
    #pragma omp single nowait
    std::cout << fibonacci(n) << std::endl;
  }
}