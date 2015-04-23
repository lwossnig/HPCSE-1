#include "mpi.h"

int main(int argc, char** argv){
	MPI_Init(&argc, &argv);
	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(rank == 0){
		MPI_Status status;
		char txt[100];
		MPI_Recv(txt, 100, MPI_CHAR, 1, 42, MPI_COMM_WORLD, &status);
		std::cout << txt << "\n";
	}
	else{
		std::string text="Hello world!";
		MPI_Send(const_cast<char*>(text.c_str()), text.size()+1, MPI_CHAR, 0, 42, MPI_COMM_WORLD); 
	}
	// e.g. Tag: MPI_ANY_TAG
	//


	MPI_Finalize();

	return 0;
}

