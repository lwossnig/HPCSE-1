#include <iostream>
#include <vector>
#include <numeric>
#include <cassert>
#include <mpi.h>

int main(int argc, char* argv[])
{
	// vector size
	const int N = 1600000;
	int num_processes, rank;
	
	// SET UP COMMUNICATION
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// DETERMINE num_processes AND OUR rank

	// initialize local parts of the vectors and do the sum z = x + y
	assert( N % num_processes == 0 );
	int nlocal = N / num_processes;
	std::vector<float> x(nlocal,-1.2), y(nlocal, 3.4), z(nlocal);
	for( int i = 0; i < nlocal; i++ )
		z[i] = x[i] + y[i];

	float * fullz = NULL;

	if ( rank == 0 )
	{
		std::vector<float> fullz(N);
		MPI_Status status;

		// COLLECT ALL PARTS INTO fullz

		for( size_t i = 0; i < nlocal; i++ )
			fullz[i] = z[i];

		for(size_t i = 1; i < num_processes ; ++i)
		{
			MPI_Recv(&fullz[i*nlocal], nlocal,MPI_FLOAT, i, 42, MPI_COMM_WORLD, &status);
		}


		//print result checksum
		std::cout << std::accumulate( fullz.begin(), fullz.end(), 0. ) << std::endl;

	} else {
		// SEND LOCAL PART z TO ROOT PROCESS
		MPI_Send(&z[0], nlocal, MPI_FLOAT, 0, 42, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}
