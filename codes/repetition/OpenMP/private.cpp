#include <iostream>
#include <omp.h>

int k = 8;
#pragma omp threadprivate (k)
int main()
{
	int i = 9;
	int n = 100;
	int sum = 0;
#pragma omp parallel firstprivate(i) 
	{
	std::cout << omp_get_num_threads() << std::endl;
	i += 1;
	std::cout << "My thread num = " << omp_get_thread_num() << std::endl;
	std::cout << "In thread i:" << i << std::endl;
	std::cout << "In thread n:" << n << std::endl;
	sum += i;
	std::cout << "In thread k:" << k << std::endl;
	}
#pragma omp single
	std::cout << "i after : " << i << std::endl;
	std::cout << "sum after : " << sum << std::endl;

	return 0;
}
