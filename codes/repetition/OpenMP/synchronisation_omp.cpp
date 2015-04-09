#include <iostream>
#include <omp.h>

int main()
{
#pragma omp parallel
	{
#pragma omp critical (output1)
		std::cout << "I am thread " << omp_get_thread_num() << " of " << omp_get_num_threads() - 1 << " threads." << std::endl;
#pragma omp barrier
#pragma omp single
		std::cout << "Blablabla.... I'm ("<< omp_get_thread_num() <<") waiting to go on!" << std::endl; /// no two threads will simultaneously be in a critical section with the same name (outputx)
	}
}
