/// Code based on HPCSE14, Matthias Troyer
/// Lecture and Exercises, see www.cse-lab.ethz.ch
/// (C) Leonard Wossnig for this implementation
#include <iostream>
#include <thread>
#include <mutex>
#include <cassert>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

#include "barrier.hpp"
#include "timer.cpp"


namespace diffusion
{

typedef double value_type;
typedef std::vector<value_type> container_type;

class Diffusion
{
	public:
		/// constructor
		Diffusion(const value_type D,
			  const value_type L, 
			  const value_type N, 
			  const	value_type dt,
			  const size_t nthreads)
		: D_(D), L_(L), N_(N), Ntot_(N_*N_), dt_(dt), nthreads_(nthreads), b_(nthreads)
		{
			dr_ = L_ / (N_ - 1); /// real space grid
			fac_ = D_ * dt_ / (dr_ * dr_); /// stencil factor

			rho.resize(Ntot_,0); /// initialize the lattice
			rho_tmp.resize(Ntot_,0);

			initial_density(); /// intializing starting density and dirichlet bc
		};

		void advance(const size_t start, const size_t stop)
		{ /// using Dirichlet boundaries, forward euler in time and central differences in space
			for(size_t i= start; i < stop; ++i){
				for(size_t j=0; j < N_; ++j){
					rho_tmp[i*N_+j] = rho[i*N_+j]
						+ fac_ *
						(
						 (j== 0 ? 0. : rho[i*N_ + (j - 1)]) 
						 +
						 (j== N_-1 ? 0. : rho[i*N_ + (j + 1)]) 
						 +
						 (i== 0 ? 0. : rho[(i-1)*N_ + j]) 
						 +
						 (i== N_-1 ? 0. : rho[(i+1)*N_ + j])
						 -
						 4. * rho[i*N_ + j]
						);
				}
			}
			/// waits until all threads call wait();
			b_.wait();
			using std::swap;
			if(start == 0) /// only thread 0 does the swap (one thread)
				swap(rho_tmp, rho);

			b_.wait(); /// wait till swap is done and 0 calls wait (then all others may continue)
		}

		/// function to print the points to file
		void operator() ( std::string const& filename) const
		{
			std::ofstream out_file(filename, std::ios::out);
			//out_file << "x,y,z\n"; 
			for(size_t i = 0; i < N_; ++i) {
				for(size_t j = 0; j < N_; ++j)
					out_file << (i*dr_ - L_/2.) << " " << (j*dr_ - L_/2.) << " " << rho[i*N_ + j] << "\n";
			//	out_file << "\n";
			}
			out_file.close();
		}

	private:
		void initial_density() /// initialize initial density distribution
		{
			value_type bound = 1./2.;
			for(size_t i = 0; i < N_ ; i++){
				for (size_t j=0; j < N_ ; j++){
					if(std::abs(i * dr_ - L_/2.) < bound && std::abs(j * dr_ - L_/2.) < bound)
						rho[i*N_+j] = 1.; /// central parts set to 1 for |x,y| < 0.5, 0 otherwise; note: first multiply by N_ (i ) then sum j! performance
				}
			}
		}

		value_type D_, L_;
		value_type dr_, dt_, fac_;

		size_t N_, Ntot_;
		size_t nthreads_;

		container_type rho;
		container_type rho_tmp;

		barrier b_;
};
} // end namespace

using namespace diffusion;
int main(int argc, char* argv[])
{
	if (argc < 6) {
		std::cerr << "Usage: " << argv[0] << " D L N dt nthreads" << std::endl;
		return 1;
	}

    const value_type D  = std::stod(argv[1]);
    const value_type L  = std::stod(argv[2]);
    const size_t  N  = std::stoul(argv[3]);
    const value_type dt = std::stod(argv[4]);
    const size_t nthreads = std::stoul(argv[5]);

    std::vector<std::thread> threads(nthreads);
    
    Diffusion system(D, L, N, dt, nthreads);
    system("density_0.dat");

    const value_type tmax = 10000 * dt; 
    
    Timer t;
    char filename[160];
    int count=0;
    t.start();
    for(size_t n = 0; n<nthreads; n++){
	    threads[n] = std::thread([&,n] () {
			    const size_t start = n * N/nthreads;
			    const size_t stop = (n+1) * N/nthreads;
			    value_type time = 0;
			    while (time < tmax) {
			    count++;
			    system.advance(start, stop);
			    /*
			       if(count <= 1000){
			       sprintf(filename,"density_%d.dat", count); // now same as "picture#nr"
			       system(filename);
			       }
			       */
			    time += dt;
			    }
			    });
    }

    for(auto& t :threads)
	    t.join();

    t.stop();
    
    std::cout << "Timing : " << N << " " << t.duration() << std::endl;
    
    system("density_serial.dat");
    
    return 0;
}

