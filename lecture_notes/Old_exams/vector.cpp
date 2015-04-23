#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

int main()
{
	std::mt19937 rng;
	rng.seed(time(NULL));
	std::uniform_real_distribution<double> ureal_d(0,100);
	std::vector<double> i(10);
	auto rand = std::bind(ureal_d, rng);
	//std::vector<double> j[10];
	//std::generate(j.begin(), j.end(), rng);
	std::generate(i.begin(), i.end(), rand);
	//for(auto v: j)
	//	std::cout<< v<<std::endl;
	for(auto v: i)
		std::cout<< v<<std::endl;

	return 0;
}

