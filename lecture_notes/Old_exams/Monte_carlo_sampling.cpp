/// monte carlo sampling of the integral:
/// f(x1, ..., x4) = Product(i=1;N) sin(xi + pi/2 * i)
/// in the box xi element of [-2, 3/2)

#include <iostream>
#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <omp.h>

#define DIM 4

inline double func(std::array<double,DIM> const& x)
{
	double ret=1.;
	for (std::size_t i=0; i<DIM; ++i)
		ret *= std::sin(x[i] + M_PI/2. * (i+1));
	return ret;
}

int main(int argc, char* argv[])
{
	const std::size_t M = 10000000;
	const double a=-2., b=1.5;
	long double avgx2;

	const double volume = std::pow((b-a), DIM);

	std::mt19937 eng(11);
	std::uniform_real_distribution<double> dist(a,b);
	auto rng = std::bind(dist, std::ref(eng));

	long double sum=0.;


#pragma omp parallel for reduction(+:sum,avgx2) schedule(static)
	for (std::size_t m=0; m<M; ++m)
	{
		
		std::array<double, DIM> x;
		std::generate(x.begin(), x.end(), rng);
		long double val = func(x);
		avgx2 += val*val;
		sum += val;
	}

	double mean = sum / double(M);
	avgx2 /= double(M);
	
	long double error =(avgx2 - mean*mean)/(double(M-1.)) ;
	/// here : variance = N/(N-1) * ( average(X^2) - average(X)^2 )
	
	double integral = volume * mean;
	std::cout << "Result is " << integral << " +/- " << std::sqrt(error) * volume  << std::endl;
}
