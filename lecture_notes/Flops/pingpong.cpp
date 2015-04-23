// pingpong.cpp
// ping pong send a buffer between two mpi nodes and report the times

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "mpi.h"

int main(int argc, char **argv){
    int rank, processes;

    // init MPI and find our rank and the total amount of processes
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    // assuming processes on the same node are numbered consecutively
    // (which might not be the case) this should hopefully give us a
    // process on a different node. So we see actual network bandwidth
    // as opposed to intranode bandwidth which might be optimized for.
    int client = processes-1;

    if(rank == 0) {
        std::cout << '#' << std::setw(20) << "size [bytes]" << std::setw(20) << "duration [ms]" << std::endl;
    }

    for(int s = 16;s<=128*1024*1024;s*=2)
    {
        // make sure all nodes have caught up
        MPI_Barrier(MPI_COMM_WORLD);

        std::vector<char> buffer(s);
        size_t byte_size = sizeof(char)*buffer.size();

        int iterations = std::min(1024*1024*1024/s, 1024);

        double start = MPI_Wtime();
        // ping pong the buffer iterations times.
        for(int j = 0;j<iterations;++j) {
            MPI_Status status;
            if(rank == client) {
                MPI_Recv(&buffer[0], byte_size, MPI_CHAR, 0, j, MPI_COMM_WORLD, &status);
                MPI_Send(&buffer[0], byte_size, MPI_CHAR, 0, j, MPI_COMM_WORLD);
            } else if(rank == 0) {
                MPI_Send(&buffer[0], byte_size, MPI_CHAR, client, j, MPI_COMM_WORLD);
                MPI_Recv(&buffer[0], byte_size, MPI_CHAR, client, j, MPI_COMM_WORLD, &status);
            }
        }
        double end = MPI_Wtime();

        // estimate the duration of one transfer by averaging
        double duration = (end-start)/(2*iterations)*1e3;

        if(rank == 0) {
            std::cout << std::setw(20) << byte_size << std::setw(20) << duration << std::endl;
        }
    }

    // shutdown MPI
    MPI_Finalize();

    return 0;
}
