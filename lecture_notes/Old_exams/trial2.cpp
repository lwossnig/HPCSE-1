#include <iostream>
#include <omp.h>

int x;
int main()
{
	x=0;
#pragma omp parallel shared(x)
	{
#pragma omp critical
		x += 1;
	}
	std::cout << x << std::endl;
	return 0;
}
