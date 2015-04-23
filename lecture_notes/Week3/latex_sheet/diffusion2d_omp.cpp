// CODE BASED ON THE LECTURE (solution) CODE FROM HPCSE-I, ETHZ
// CHANGES ARE ACCORDINGLY JUST DONE BY THE PARALLELIZATION 
// OF USING OPEN MP
// Leonard Wossnig, 14.10.14
// compile and link using:
// g++ -fopenmp <codename>.cXX
// using #pragram omp parallel 


#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <omp.h>
#include <thread>
#include <cstdlib>
#include "timer.cpp"


typedef double value_type;
typedef std::size_t size_type;

class Diffusion2D {

public:
    Diffusion2D(const value_type D, //Diffusion constant
                const value_type L, // length of interval
                const size_type N,  //steps of grid (1d)
                const value_type dt) //timestep
    : D_(D), L_(L), N_(N), Ntot(N_*N_), dt_(dt)
    {
        /// real space grid spacing
        dr_ = L_ / (N_ - 1);
        
        /// stencil factor
        fac_ = dt_ * D_ / (dr_ * dr_);
        
        rho_.resize(Ntot, 0.);
        rho_tmp.resize(Ntot, 0.);
        
        initialize_density();
	
    }
    
    void advance()
    {
        /// Dirichlet boundaries; central differences in space, forward Euler
        /// in time
#pragma omp parallel for collapse(2)	    
        for(size_type i = 0; i < N_; ++i) {
            for(size_type j = 0; j < N_; ++j) {
                rho_tmp[i*N_ + j] = rho_[i*N_ + j] +
                fac_
                *
                (
                 (j == N_-1 ? 0. : rho_[i*N_ + (j+1)])
                 +
                 (j == 0    ? 0. : rho_[i*N_ + (j-1)])
                 +
                 (i == N_-1 ? 0. : rho_[(i+1)*N_ + j])
                 +
                 (i == 0    ? 0. : rho_[(i-1)*N_ + j])
                 -
                 4.*rho_[i*N_ + j]
                 );
            }
	} 
        /// use swap instead of rho_=rho_tmp. this is much more efficient, because it does not have to copy
        /// element by element.
            
#pragma omp barrier
#pragma omp master
	using std::swap;
        swap(rho_tmp, rho_);
//#pragma omp flush 
    }
    
    void write_density(std::string const& filename) const
    {
        std::ofstream out_file(filename, std::ios::out);
        
        for(size_type i = 0; i < N_; ++i) {
            for(size_type j = 0; j < N_; ++j)
                out_file << (i*dr_ - L_/2.) << '\t' << (j*dr_ - L_/2.) << '\t' << rho_[i*N_ + j] << "\n";
            out_file << "\n";
        }
        out_file.close();
    }

    
    value_type mc_rho()
    {
	    size_type mc_steps2_=0;
    	    value_type mc_result2_=0, mc_sum2_=0;

	    while(mc_steps2_ < 10000)
	    {
		    mc_sum2_+=rho_[(rand()%N_)*N_ + (rand()%N_)];
		    mc_steps2_++;
	    }
	    mc_result2_ = L_*mc_sum2_/10000.;
	    return (mc_result2_);
    }


    value_type mc_cloud()
    {	    
	    size_type mc_steps_=0;
    	    value_type mc_result_=0, mc_sum_=0;

	    size_type j=0;
	    size_type i=0;
	    while(mc_steps_ < 10000)
	    {
		    j = (rand()%N_); //really bad rng but quick implementation for now
		    i = (rand()%N_); // -"-
		    mc_sum_ += rho_[i*N_+j]*((i*dr_ - L_/2.)*(i*dr_ - L_/2.) + (j*dr_ - L_/2.)*(j*dr_ - L_/2.)) ;
		    mc_steps_++;
	    }
	    mc_result_ = L_*mc_sum_/10000.;
	    return (mc_result_);
    }
    
private:
    
    void initialize_density()
    {
        /// initialize rho(x,y,t=0)
        value_type bound = 1/2.;
        
        for (size_type i = 0; i < N_; ++i) {
            for (size_type j = 0; j < N_; ++j) {
                if (std::abs(i*dr_ - L_/2.) < bound && std::abs(j*dr_ - L_/2.) < bound) {
                    rho_[i*N_ + j] = 1;
                } else {
                    rho_[i*N_ + j] = 0;
                }
            }
        }
	std::cout << "Initialized density, now running programm with " << std::thread::hardware_concurrency() <<" concurrent supported threads in parallel.\n";
#pragma omp parallel
#pragma omp critical (output) 
	std::cout << "Thread " << omp_get_thread_num() << " of " << omp_get_num_threads() << " threads." << std::endl;

    }
    
    value_type D_, L_;
    size_type N_, Ntot;
    
    value_type dr_, dt_, fac_;
// MC variables
    
    std::vector<value_type> rho_, rho_tmp;
};


int main(int argc, char* argv[])
{
    std::ofstream out("rho_evolution.txt");
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " D L N dt" << std::endl;
        return 1;
    }
    
    const value_type D  = std::stod(argv[1]);
    const value_type L  = std::stod(argv[2]);
    const size_type  N  = std::stoul(argv[3]);
    const value_type dt = std::stod(argv[4]);

    
    Diffusion2D system(D, L, N, dt); 
    system.write_density("density.0.dat");

    srand(time(NULL));

    const value_type tmax = 10000 * dt;
    value_type time = 0;
    int w_count= 0; 
    Timer t;
    t.start(); 

    while(time < tmax){
        system.advance();
	w_count++;
	if(w_count%100==0){
	  out << time << "    " << system.mc_rho() << "   " << system.mc_cloud() <<std::endl;
	}
	time += dt;
    }

    t.stop();
    
    std::cout << "Timing : " << N << " " << t.duration() << std::endl;
    
    system.write_density("density_serial.dat");
   
    out.close();
    return 0;
}
