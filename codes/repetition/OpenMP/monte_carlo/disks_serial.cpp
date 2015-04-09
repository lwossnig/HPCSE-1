
#include <iostream>
#include <vector>
#include <string>
#include <omp.h>
#include "disks_mc.hpp"
#include "vector_stats.hpp"
#include "timer.hpp"


int main () {
    std::size_t Nx = 14;
    std::size_t Ny = 16;
    double L = 1.;
    std::size_t nu = 5;
    double d0 = (1. - (std::pow(2., nu-8.))) * L / Nx;



    std::size_t M = 512;
    double dA2 = 2.*(L/2.)*(L/2.) / M;
    
    std::size_t nequi = 64;
    std::size_t nmeas = 64;
    
    timer t;
    t.start();
   
    int nthreads;
#pragma omp parallel 
#pragma omp master
    nthreads = omp_get_num_threads();
    std::vector<VectorStats> thread_accumulators(nthreads, VectorStats(M));

#pragma omp parallel num_threads(2) 
    {
	    int thread_num =omp_get_thread_num();
	    DisksMC system(Nx, Ny, L, d0, thread_num);


	    /// Equilibration
	    for (std::size_t n=0; n<nequi; ++n)
		    system.sweep();
    
	    /// Measurements
#pragma omp for
	    for (std::size_t n=0; n<nmeas; ++n) {
		    system.sweep();
		    thread_accumulators[thread_num] << system.distance_histogram(M, dA2);
	    }
    }
    
    t.stop();
    std::cout << "# Time elapsed : " << t.get_timing() << " nthreads: " << nthreads << std::endl;

}
