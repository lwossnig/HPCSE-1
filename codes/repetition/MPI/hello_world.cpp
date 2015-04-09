#include <iostream>
#include <mpi.h>
#include <cmath>

int main(int argc, char* argv[]){

	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int flag=0;MPI_Initialized(&flag);
	if (flag)
		std::cerr << "MPI INITIALIZED FOR PROCESS " << rank <<"!" << std::endl;


	if (rank==0)  { /// Master process
		MPI_Status status;
		char txt[100];
		MPI_Recv(&txt[0], 100, MPI_CHAR, MPI_ANY_SOURCE, 42, MPI_COMM_WORLD, &status);
		std::cout << txt << std::endl;
	}
	else  { /// Worker processes
		std::string text = "Hello world!";
		MPI_Send(const_cast<char*>(text.c_str()), text.size()+1, MPI_CHAR, 0, 42, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	int flag2 = 0;MPI_Finalized(& flag2);
	if(flag2)
		std::cout << "MPI FINALIZED!\n";

	return 0;
}
