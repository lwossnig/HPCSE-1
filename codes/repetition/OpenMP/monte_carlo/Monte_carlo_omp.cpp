#include <iostream>
#include <omp.h>
#include <cassert>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <random>
#include "timer.cpp"
#define PARALLEL

namespace montecarlo
{

typedef double value_type;
typedef std::vector<value_type> container_type;

class MC 
{
	public:
		/// constructor
		MC(
			  value_type x ,
			  value_type y ,
			  value_type d ,
			  size_t steps
			  )
		: x0_(x), y0_(y) , x_(x0_), y_(y0_), d_(d), steps_(steps)
		{
			measure_.resize(steps);
		};

		void mc_steps()
		{ /// random walk from position x0_ y0_
			value_type result=0;
#ifdef PARALLEL
#pragma omp parallel for reduction(+:result)
#endif
			for (size_t i = 0; i < steps_; i++){
				std::mt19937 rng(i);
				std::uniform_real_distribution<value_type> step_distribution(0.,1.);
				x_ = x0_; y_ = y0_;
				value_type a;
				while ( check_location( x_, y_) )
				{
					a = step_distribution(rng);
					x_ += d_ * std::cos(2*M_PI*a);
					y_ += d_ * std::sin(2*M_PI*a);
				}
//				std::cout << x_ << "    " << y_ << std::endl;
				measure_[i] = evaluate(x_,y_);
				result += measure_[i];
			}
#ifdef PARALLEL
#pragma omp single
#endif
			std::cout << "The average value for rho = " << result/(double)steps_ << std::endl;

		}

		value_type evaluate(value_type x, value_type y)
		{
			return ( x < 0 ? 0 : ( x > 1 ? 1 : x ));
		}



		bool check_location (value_type x, value_type y)
		{
			return ( x < 1 && y < 1 && x > 0 && y > 0 );
		}



		/// function to print the resulting points to file
		void operator() ( std::string const& filename) const
		{
			std::ofstream out_file(filename, std::ios::out);
			for(size_t i = 0; i < steps_; i++)
				out_file << steps_<< "\t" << measure_[i] << "\n";
		      	out_file.close();
		}

	private:
		value_type x0_;
		value_type y0_;
		value_type x_;
		value_type y_;
		value_type d_;
		size_t steps_;
		std::vector<value_type> measure_;
};
} // end namespace

using namespace montecarlo;
int main(int argc, char* argv[])
{
	if(argc<2){
		std::cerr << "Usage: " << argv[0] << " steps" <<std::endl;
	}
	char filename[160];
	size_t steps = std::stoul(argv[1]);
	sprintf(filename, "measurement_%d.dat",(int)steps);
	
	Timer t;
	t.start();
	MC system(0.3, 0.4, 0.01, steps);
	system.mc_steps();
	t.stop();
	system(filename);
	std::cout << "time: " << t.duration() << std::endl;
	return 0;
}

