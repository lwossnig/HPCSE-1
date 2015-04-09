#include <iostream>
#include <cassert>
#include <new>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <complex>
//#define DEBUG
#define pi 3.14159265359

void update(std::complex<double> * , std::complex<double> * , const int , const double , const double, int, int , MPI_Status);// update using the velocity and locality vector
void write_to_file(const char * , std::complex<double> * , int N);

int main(int argc, char** argv)
{
	MPI_Status status;
	MPI_Init(&argc, &argv);
	int myrank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int N = 1000;
	double h = 1/(double)N;
	std::cout << " Starting ... " << std::endl;
	std::complex<double>* vortices = new std::complex<double>[N];
	std::complex<double>* velocities = new std::complex<double>[N];
	
	std::cout << " Initialize particle positions... " << std::endl;
	for(int i=0;i<N;i++){
		vortices[i] = {-0.5+(i+0.5)*h,0};
		velocities[i] = {0,0};
#ifdef DEBUG
		std::cout << vortices[i] << std::endl;
#endif
	}
	double t = 0; const double dt = 0.001; int no = 0, count = 0;
	char filename[160];
	/// updating the positions now
	while(t<1){
		if(myrank==0){
			if(count%100==0){
				sprintf(filename,"vortices%d.txt",no); // now same as "picture#nr"
				write_to_file(filename, vortices, N);
				no++;
			}
		}
		update(vortices, velocities, N, h, dt, myrank, size, status);
#ifdef DEBUG
		for(int i=0;i<N;i++){
			std::cout << vortices[i] << std::endl;
		}
#endif
		count++;
		std::cout << "iteration with time t = " << t << std::endl;
		t += dt;

	}
	std::cout << " System in final state. t = " << t << std::endl;
	delete [] vortices, velocities;
	MPI_Finalize();
	return 0;
}


void write_to_file(const char* filename, std::complex<double> *positions, int N){
	std::ofstream out(filename);
	for(int i=0; i < N;i++){
		out << positions[i].real() << "    " << positions[i].imag() << std::endl;
	}
	out.close ();
}

void update(std::complex<double> *positions, std::complex<double> *velocities, const int N, const double h, const double dt,int myrank, int size, MPI_Status status)
{
	std::complex<double> x(0);
	std::complex<double> ic(0,1);
	
	int tag = 42;
	int chunk_size = N/size;
	assert(N%size==0);
	
	for(int i=myrank*chunk_size; i<( (myrank+1) * chunk_size); ++i){
		velocities[i] = {0,0};
		for(int j=0; j<N; ++j){
			if(j!=i){ 
				x = {positions[j].real(),0};
				velocities[i] += (-1.) * 4. * x * h/std::sqrt(1.-4.*x*x) * ic / (2. * pi * ( positions[i] - positions[j] )); 
#ifdef DEBUG
				std::cout << velocities[i] << std::endl;
#endif
			}
		}
		positions[i] =  positions[i] + dt * std::conj(velocities[i]);
	}

	/// communication with other parts
	int sbuf_size = sizeof(std::complex<double>)*chunk_size + MPI_BSEND_OVERHEAD;
	char * buffer = new char[sbuf_size];
	MPI_Buffer_attach(buffer, sbuf_size);
	for(int com = 0; com < size; com++){
		MPI_Status status2;
		/// allocate a buffer and attach to MPI
		int from = ((myrank-com)%size);
		if(from < 0){
			from = size + from;
		}
		int to = std::abs((myrank+com)%size);
		MPI_Bsend(&positions[(myrank*chunk_size)], chunk_size, MPI_DOUBLE_COMPLEX,to,tag,MPI_COMM_WORLD);
		MPI_Recv (&positions[from*chunk_size], chunk_size, MPI_DOUBLE_COMPLEX, from, tag, MPI_COMM_WORLD, &status2);
		
		/// detach the buffer making sure all sends are done
	}
	MPI_Buffer_detach(buffer, &sbuf_size);
	delete [] buffer;
	/// communication done!
}



