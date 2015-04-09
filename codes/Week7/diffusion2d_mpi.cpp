/*(C) Leonard Wossnig - Exercise 7 - HPCSE 14 @ ETHZ */

#include <iostream>
#include <mpi.h>
#include <algorithm>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>

#define _DEBUG

#include "timer.cpp"


typedef double value_type;
typedef std::size_t size_type;

class Diffusion2D {
public:
    Diffusion2D(const value_type D,
                const value_type L,
                const size_type N,
		const size_type myrank,
		const size_type num_threads,
                const value_type dt)
    : D_(D), L_(L), N_(N), Ntot(N_*N_), dt_(dt), myrank_(myrank), num_threads_(num_threads)
    {
        /// real space grid spacing
        dr_ = L_ / (N_ - 1);
        
        /// stencil factor
        fac_ = dt_ * D_ / (dr_ * dr_);
	assert(N%num_threads_==0);
        interval_ = N_/num_threads_;
        rho_.resize(Ntot, 0.);
        rho_tmp.resize(Ntot, 0.);
        
        initialize_rho_();
#ifdef _DEBUG
	std::cout << "Process " << myrank << " started calulating the intervale [" << myrank*interval_ <<"," << (myrank+1)*interval_ << "]\n";
#endif
    }
    
    
     

    void advance()
    {
	    /// Dirichlet boundaries; central differences in space, forward Euler
	    /// in time

	    /// decompose the whole domain into N_/num_threads layers
	    /// amd each thread working on the range of [N_/num_threads *myrank, N/num_threads * (myrank+1)] 

	    for(size_type i = myrank_*interval_ ; i < (myrank_+1)*interval_; ++i) {
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
	    using std::swap;
	    swap(rho_tmp, rho_);

	    /// Now exchanging the rows of the upper and lower phantom

	    if(myrank_%2==0 && myrank_>0 && myrank_ < num_threads_ - 1){
		    int left = myrank_ - 1;
		    int right = myrank_ + 1;
#ifdef _DEBUG
		    std::cout << "Process " << myrank_ << " sending to " << left << std::endl;
#endif
		    MPI_Send(&rho_[myrank_*interval_],interval_,MPI_DOUBLE,left,0,MPI_COMM_WORLD);
		    MPI_Recv(&rho_[left*interval_],interval_,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&status);
		    MPI_Send(&rho_[myrank_*interval_],interval_,MPI_DOUBLE,right,0,MPI_COMM_WORLD);
		    MPI_Recv(&rho_[right*interval_],interval_,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&status);
	    }

	    /// size = max_rank+1 --> to exclude rank N from interaction we need to exclude size - 2 !
	    if(myrank_%2!=0 && myrank_>0 && myrank_ < num_threads_ - 2)
	    {
		    int left = myrank_ - 1;
		    int right = myrank_ + 1;
#ifdef _DEBUG
		    std::cout << "Process " << myrank_ << " receiving from " << right << std::endl;
#endif
		    MPI_Recv(&rho_[right*interval_],interval_,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&status);
		    MPI_Send(&rho_[myrank_*interval_],interval_,MPI_DOUBLE,right,0,MPI_COMM_WORLD);
		    MPI_Recv(&rho_[left*interval_],interval_,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&status);
		    MPI_Send(&rho_[myrank_*interval_],interval_,MPI_DOUBLE,left,0,MPI_COMM_WORLD);
	    }

	    /// send recv for last thread (boundary exception)
	    if(myrank_ == num_threads_ - 1){
		    int left = myrank_ - 1;
		    if(myrank_%2==0){
			    std::cout << "Process " << myrank_ << " sending to " << left << std::endl;
			    MPI_Send(&rho_[myrank_*interval_],interval_,MPI_DOUBLE,left,0,MPI_COMM_WORLD);
			    MPI_Recv(&rho_[left*interval_],interval_,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&status);
		    }
		    else{
			    std::cout << "Process " << myrank_ << " receiving from " << left << std::endl;
			    MPI_Recv(&rho_[left*interval_],interval_,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&status);
			    MPI_Send(&rho_[myrank_*interval_],interval_,MPI_DOUBLE,left,0,MPI_COMM_WORLD);
		    }
	    }

	    /// send recv for first thread (boundary exception)
	    if(myrank_ == 0){
		    int right = myrank_ + 1;
		    std::cout << "Process " << myrank_ << " sending to " << right << std::endl;
		    MPI_Send(&rho_[myrank_*interval_],interval_,MPI_DOUBLE,right,0,MPI_COMM_WORLD);
		    MPI_Recv(&rho_[right*interval_],interval_,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&status);
	    }

	    
    }

    /// function to summarize all results in the master process from all slaves! 
    void final_rho(){
	    if(myrank_==0){
		    for(size_type i =1; i<num_threads_;i++){
			    MPI_Recv(&rho_[i*interval_],interval_,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
		    }
	    }
	    else{
		    /// here i chose synchronous send to assure the safe transmission to the master process (speed doesn't matter here anymore!)
		    MPI_Ssend(&rho_[myrank_*interval_],interval_,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	    }
    }




    void write_rho_(std::string const& filename) const
    {
        std::ofstream out_file(filename, std::ios::out);
	out_file << "x,y,z\n"; 
        for(size_type i = 0; i < N_; ++i) {
            for(size_type j = 0; j < N_; ++j)
                out_file << (i*dr_ - L_/2.) << "," << (j*dr_ - L_/2.) << "," << rho_[i*N_ + j] << "\n";
            out_file << "\n";
        }
        out_file.close();
    }
    
private:
    
    void initialize_rho_()
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
    }
    
    value_type D_, L_;
    size_type N_, Ntot;
    size_type myrank_;
    value_type interval_;
    size_type num_threads_;
    value_type dr_, dt_, fac_;
    MPI_Status status; 
    std::vector<value_type> rho_, rho_tmp;
};


int main(int argc, char** argv)
{
    
    MPI_Init(&argc, &argv);
    int myrank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifdef _DEBUG
    std::cout << "Process " << myrank << " of " << size << " initialized\n";
#endif


    if (argc < 5) {
	    std::cerr << "Usage: " << argv[0] << " D L N dt" << std::endl;
	    return 1;
    }
    const value_type D  = std::stod(argv[1]);
    const value_type L  = std::stod(argv[2]);
    const size_type  N  = std::stoul(argv[3]);
    const value_type dt = std::stod(argv[4]);
    
    /* e.g.:  
    value_type D  = 1;
    value_type L  = 1;
    size_type N = 128;
    value_type dt = 0.00001;
    */
    
    /// on each thread initialize whole system!
    Diffusion2D system(D, L, N, myrank, size, dt);
    
    /// only master thread!
    if(myrank==0){
	    system.write_rho_("density.0.csv");
    }

    const value_type tmax = 10000 * dt; 
    value_type time = 0;
    
    Timer t;
    
    t.start();
    while (time < tmax) {
        system.advance();
#ifdef _DEBUG
	std::cout << "Time: " << time << std::endl;
#endif
        time += dt;
	system.final_rho();
    }
    t.stop();
    
    std::cout << "Timing : " << N << " " << t.duration() << std::endl;
    
    /// Master thread should receive all the arrays from the slave processes here!
    /// should only be done by master!
    if(myrank==0){
	    system.write_rho_("density_serial.csv");
    }

    MPI_Finalize();
    return 0;
}
