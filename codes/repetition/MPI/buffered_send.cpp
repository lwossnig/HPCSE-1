#include <iostream>
#include <vector>
#include <mpi.h>

typedef std::vector<double> container;
int main (int argc, char* argv[]) {
	MPI_Status status;
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	size_t n =  100;

	container ds(n,1) ;
	container dr(n,0) ;
	int tag=99;
	if(rank == 0){
		std::cout << "Process " << rank << " dr[0] = " << dr[0] << std::endl;
	}

	/// allocating buffer and attach it to MPI
	int buffer_size = sizeof(double) * n + MPI_BSEND_OVERHEAD;
	char* buffer = new char[buffer_size];
	MPI_Buffer_attach(buffer, buffer_size);


	if (rank==0)  {
		MPI_Bsend (&ds[0], n, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
		MPI_Recv  (&dr[0], n, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD, &status);
	}

	else	      {
		MPI_Bsend (&ds[0], n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
		MPI_Recv  (&dr[0], n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
	}


	/// detach the buffer
	MPI_Buffer_detach(buffer, &buffer_size);
	delete[] buffer;

	if(rank == 0){
		std::cout << "Process " << rank << " dr[0] = " << dr[0] << std::endl;
	}

	MPI_Finalize();

	return 0;
}



