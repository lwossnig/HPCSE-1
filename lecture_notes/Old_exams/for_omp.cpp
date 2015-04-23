#include<iostream>
#include <omp.h>
#include <vector>
#include <random>
#include "timer.hpp"
  
int main()
{
	int N=40;
	std::vector<int> result1(40,0);
	std::vector<int> result2(40,2);
	std::vector<int> B(40,1);
	std::vector<int> A(40,1);
	std::vector<int> idx(40,1);
	std::mt19937 mt;
	std::uniform_int_distribution<int> uint_d(0,40);
	mt.seed(42);

	for (int i=0; i<40; i++){
		idx[i] = uint_d(mt);
	}

	for (int i=0; i<40; i++){
		result1[i] = i;
		B[i] = i*i;
		std::cout<< result1[i] << "\t";
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
#ifdef PARALLEL
#pragma omp parallel
	std::cout <<"Here is thread: " << omp_get_thread_num()+1 << " of " << omp_get_num_threads() << std::endl;
#pragma omp parallel for
#endif
	/******* HERE TEST AREA STARTS *******/
	for (int i=1; i<N-1; i++){
		result1[i] = 2*A[i];
		result2[i] = result1[i-1]+B[i];
	}

	for (int i=0; i<40; i++)
		std::cout<< result1[i] << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (int i=0; i<40; i++)
		std::cout<< result2[i] << "\t";
	std::cout << std::endl;
	std::cout << std::endl;

	return 0;
}


