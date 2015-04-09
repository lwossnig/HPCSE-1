/// Code based on HPCSE14, Matthias Troyer
/// Lecture and Exercises, see www.cse-lab.ethz.ch
/// (C) Leonard Wossnig for this implementation
#include <iostream>
#include <omp.h>
#include <cassert>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

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
			  const	value_type dt)
		: D_(D), L_(L), N_(N), Ntot_(N_*N_), dt_(dt)
		{
			dr_ = L_ / (N_ - 1); /// real space grid
			fac_ = D_ * dt_ / (dr_ * dr_); /// stencil factor

			rho.resize(Ntot_,0); /// initialize the lattice
			rho_tmp.resize(Ntot_,0);

			initial_density(); /// intializing starting density and dirichlet bc
		};

		void advance()
		{ /// using Dirichlet boundaries, forward euler in time and central differences in space
#pragma omp parallel for collapse(2)
			for(size_t i= 0; i < N_; ++i){
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
			using std::swap;
			swap(rho_tmp, rho);
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

		container_type rho;
		container_type rho_tmp;
};
} // end namespace

using namespace diffusion;
int main(int argc, char* argv[])
{
	if (argc < 5) {
		std::cerr << "Usage: " << argv[0] << " D L N dt " << std::endl;
		return 1;
	}
#pragma omp parallel
	{
#pragma omp master
		std::cout << "Running with " << omp_get_num_threads() <<" threads" << std::endl;
	}

    const value_type D  = std::stod(argv[1]);
    const value_type L  = std::stod(argv[2]);
    const size_t  N  = std::stoul(argv[3]);
    const value_type dt = std::stod(argv[4]);

    Diffusion system(D, L, N, dt);
    system("density_0.dat");

    const value_type tmax = 10000 * dt; 
    
    Timer t;
    char filename[160];
    int count=0;
    value_type time = 0;
    t.start();
    while (time < tmax) {
	    count++;
	    system.advance();
	    /*
	       if(count <= 1000){
	       sprintf(filename,"density_%d.dat", count); // now same as "picture#nr"
	       system(filename);
	       }
	       */
	    time += dt;
    }
    t.stop();
    
    std::cout << "Timing : " << N << " " << t.duration() << std::endl;
    
    system("density_serial.dat");
    
    return 0;
}

