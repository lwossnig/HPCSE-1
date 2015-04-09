
#include <omp.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include "timer.cpp"
#define _PARALLEL
#define pi 3.14159265359

typedef double value_type;
typedef std::size_t size_type;


class Diffusion2D {

/*********************************************************/
/*********************************************************/
public:
    Diffusion2D(const value_type D,
                const value_type L)
    : D_(D), L_(L)
    {}
   // D unneccesary here but left in for trial of a D-dependend trial of random walks...
   
/*********************************************************/
    static value_type rho_function(value_type x, value_type y)
    {
	    return x;
    }

/*********************************************************/
    value_type advance(value_type x, value_type y, int Nin)
    {

	  
	  const value_type d = 0.01;
	  int N = Nin;
	  std::vector<int> seeds(N,0);
	  value_type sum = 0.;
	  std::random_device rd;
	  std::mt19937 mt1(rd());
	  std::uniform_int_distribution<int> uint_d(0,N);
	  for(size_t k=0; k<N; k++)
	  {
		  seeds[k] = uint_d(mt1);
	  }



#ifdef _PARALLEL
omp_set_num_threads(4);
std::cout << omp_get_num_threads() << std::endl;
#pragma omp parallel for reduction(+:sum)
#endif
	  for( size_type k=0; k<N ; k++ )
	  {  
		  value_type delta_x = x, delta_y = y;
		  std::mt19937 mt2(seeds[k]);
		  std::uniform_real_distribution<value_type> ureal_d(0.,1.);
		  // iteration random walk till border crossed
#ifndef NDEBUG
#pragma omp critical (output)
		  std::cout << "Thread " << omp_get_thread_num() << " going into while loop" << std::endl;
#endif
		  while(true)
		  {
			  delta_x += d*std::cos(2 * pi * ureal_d(mt2) );
			  delta_y += d*std::sin(2 * pi * ureal_d(mt2) );
			  if( (delta_x >= L_) || (delta_y >= L_) || (delta_x <= 0.) || (delta_y <= 0.) )  
			  {
				  break;
			  }
		  }
  
		  // if border crossed --> add rho(x_c,y_c) to sum
		  // incase i want to adjust the exact position determination i can adjust different for crossing of x and y boder
		  // but in this case it should average equally if i just take the position at the end of the crossing step
		  if( (delta_x >= L_/2.) || (delta_y >= L_/2.) || (delta_y <= 0.) )
		  {
			  sum += rho_function(delta_x, delta_y) ;
		  }
		  else//if( (delta_x <= 0.) || (delta_y <= 0.) )
		  {
			  sum += 0; 
		  }
	  }
	  return (sum/N);
    }
        

/*********************************************************/

    void write_density(value_type sum, std::string const& filename) const
    {
        std::ofstream out_file(filename, std::ios::out);
        out_file << sum << std::endl; 
        out_file.close();
    }
    

/*********************************************************/
/*********************************************************/

private:
    value_type D_, L_;
/*********************************************************/
};


int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " D L N" << std::endl;
        return 1;
    }
    
    const value_type D  = std::stod(argv[1]);
    const value_type L  = std::stod(argv[2]);
    int		     Nin  = std::stod(argv[3]);


    
    Diffusion2D system(D, L);
    value_type sum;
    Timer t;
    
    t.start();
    sum = system.advance(0.3,0.4, Nin);
    system.write_density(sum ,"density_0.3_0.4.dat");
    t.stop();
    
    //std::cout << "Timing : " << t.duration() << std::endl;
    std::cout << sum << std::endl;
    
    return 0;
}
