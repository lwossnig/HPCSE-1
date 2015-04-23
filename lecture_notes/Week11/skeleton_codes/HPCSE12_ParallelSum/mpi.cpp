// Skeleton code for HPCSE Exam, 18.12.2012
// Profs. P. Koumoutsakos and M. Troyer
// Question 5d)

#include <vector>
#include <numeric>
#include <cassert>
#include <mpi.h>
#include <iostream>

int main( int argc, char** argv )
{
    // vector size
    const int N = 1600000;
    int num_processes, rank;
    
    // SET UP COMMUNICATION.
    // DETERMINE num_processes AND OUR rank
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    
    // initialize local parts of the vectors and do the sum z = x + y
    assert( N % num_processes == 0 );
    int nlocal = N / num_processes;
    std::vector<float> x(nlocal,-1.2), y(nlocal,3.4), z(nlocal);
    for( int i = 0; i < nlocal; i++ )
        z[i] = x[i] + y[i];
    
    if( rank == 0 )
    {
        std::vector<float> fullz(N);
	MPI_Status status;        
	for(int j = 0; j<nlocal; j++)
		fullz[j]=z[j];
        
        // COLLECT ALL PARTS INTO fullz
	for(int i =1; i< num_processes; ++i)
		MPI_Recv(&fullz[i*nlocal], nlocal, MPI_FLOAT, i, 42, MPI_COMM_WORLD, &status);
        
        // print result checksum
        std::cout << std::accumulate( fullz.begin(), fullz.end(), 0. ) << std::endl;
    }
    else
    {
        // SEND LOCAL PART z TO ROOT PROCESS
	MPI_Send(&z[0], nlocal, MPI_FLOAT, 0, 42, MPI_COMM_WORLD);
    }
    
    // CLEAN UP
    MPI_Finalize();
    return 0;
}
