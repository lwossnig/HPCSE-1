#include <iostream>
#include <thread>
#include <vector>
#include <numeric>
typedef std::vector<float> container;
int main(int argc, char * argv[])
{
	// vector size
	const int N = 1600000;

	// initialize vectors
	container x(N,-1.2), y(N,3.4), z(N);
	std::vector<std::thread> threads(4);

	// DO THE SUM z = x + y using 4 threads
	std::cout << std::endl;	
	const int workload = N/4;
	for(size_t i=0; i<4; ++i)
		threads[i] = std::thread([&,i](){ 
				std::cout << "Thread " << i << " now working\n";
				for(size_t k = i * workload; k<(i+1)*workload; ++k)
				{
				      z[k] = x[k] + y[k];     
				}
				});

	for(size_t i=0; i<4; ++i)
		threads[i].join();

	// print result checksum
	std::cout<<z[1] << std::endl;
	std::cout << std::accumulate(z.begin(), z.end(), 0.) << std::endl;

//	return 0;
}

