#include <iostream>
#include <mpi.h>
#include <vector>
typedef std::vector<double> container;

int main  (int argc, char* argv[])  {

	MPI_Init (&argc, &argv);
	int rank, size;
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	MPI_Status status[2];
	MPI_Request req[2];
	size_t n = 100;

	int tag = 99;
	int left = rank-1;
	int right= rank+1;
	container in(n, 0);
	container out(n, 0);

	for (size_t i = 0; i < n; ++i){
		if  (left >=0) {
			MPI_Irecv(&in[0], n, MPI_DOUBLE, left, tag, MPI_COMM_WORLD, &req[0]);
		} else	{
			req[0] = MPI_REQUEST_NULL;
		}
		
		if  (right < size) {
			MPI_Isend(&out[0], n, MPI_DOUBLE, right, tag, MPI_COMM_WORLD, &req[1]);
		} else	{
			req[1] = MPI_REQUEST_NULL;
		}
		MPI_Waitall(2, req, status);
		in[rank] = rank;
		using std::swap;
		swap(in, out);
		
	}

	if  (rank== size-1) {
		for (size_t j=0; j < n; ++j)
			std::cout << in[j] << " , ";
		std::cout << std::endl;
	}

	MPI_Finalize();

	return 0;
}
		

