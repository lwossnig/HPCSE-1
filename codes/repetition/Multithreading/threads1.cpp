#include <iostream>
#include <vector>
#include <numeric>
#include <random>
#include <thread>

void sumterms(double& sum, std::size_t i, std::size_t j)
{
	std::cout << "thread[" << i << "] working now!" << std::endl;
	sum = 0.0;
	for(std::size_t t=i; t<j;++t){
		sum += (1.0 - 2* (t % 2)) / (2 * t + 1);
	}
}

int main()
{
	std::size_t const nthreads = std::max(1u, std::thread::hardware_concurrency());
	std::vector<double> sum(nthreads,0);
	std::vector<std::thread> threads(nthreads);
	for (unsigned i=0; i< threads.size(); i++){
		threads[i] = std::thread(
				sumterms, std::ref(sum[i]), 
				i*1000, (i+1)*1000);
	}
	for(auto& t : threads)
		t.join();
	
	double result = 0, product = 0;
	result = std::accumulate(sum.begin(), sum.end(), 0.);
	product = std::accumulate(sum.begin(), sum.end(), 1., [&](double x, double y){return x*y;});
	std::cout << result << std::endl;
	std::cout << product << std::endl;

}

