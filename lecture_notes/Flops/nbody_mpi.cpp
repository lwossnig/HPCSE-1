// nbody.cpp
// naive gravitational nbody simulation using MPI

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "mpi.h"

int main(int argc, char **argv){
    const int N = 4*1024;
    float dt = 0.001f;

    // this bool controls whether blocking or nonblocking recv is used
    // with blocking recv we should observe less transaction/compute overlap
    // and as a result worse performance.
    bool blocking = true;

    int rank, processes;
    // init MPI and find our rank and the total amount of processes
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    // we assume N is a multiple of the number of processes
    int Nlocal = N/processes;

    // position and velocity
    std::vector<float> position(3*Nlocal, rank);
    std::vector<float> velocity(3*Nlocal, rank);

    // communication buffers
    std::vector<float> otherposition(3*Nlocal, rank);
    std::vector<float> recvbuf(3*Nlocal, rank);

    // init rand differently on every node, otherwise they will all get
    // the same positions
    std::srand(rank);

    // init state (velocities are 0 per default initialisation)
    for(int i = 0;i<Nlocal;++i) {
        position[3*i+0] = 10.0f*std::rand()/RAND_MAX-5.0f;
        position[3*i+1] = 10.0f*std::rand()/RAND_MAX-5.0f;
        position[3*i+2] = 10.0f*std::rand()/RAND_MAX-5.0f;
    }

    double t_compute = 0;
    double start = MPI_Wtime();
    for(int k = 0;k<100;k++) {
        // start with own positions
        otherposition = position;
        for(int i = 1;i<processes;++i) {
            MPI_Request requests[2];
            MPI_Status statuses[2];

            // start send and receive
            int tag = i+processes*k;
            // note that using a blocking send here would result in a deadlock
            MPI_Isend(&position[0], position.size(), MPI_FLOAT, (rank+processes-i)%processes, tag, MPI_COMM_WORLD, &requests[0]);
            if(!blocking) {
                MPI_Irecv(&recvbuf[0], recvbuf.size(), MPI_FLOAT, (rank+i)%processes, tag, MPI_COMM_WORLD, &requests[1]);
            } else {
                MPI_Recv(&recvbuf[0], recvbuf.size(), MPI_FLOAT, (rank+i)%processes, tag, MPI_COMM_WORLD, &statuses[1]);
            }

            // do calculations while comm is transferring data
            double start_compute = MPI_Wtime();
            for(int i = 0;i<Nlocal;++i) {
                for(int j = 0;j<Nlocal;++j) {
                    float diffx = otherposition[3*i+0]-position[3*j+0];
                    float diffy = otherposition[3*i+1]-position[3*j+1];
                    float diffz = otherposition[3*i+2]-position[3*j+2];
                    float dist = std::sqrt(diffx*diffx + diffy*diffy + diffz*diffz)+0.001;
                    velocity[3*i+0] += dt*diffx/(dist*dist*dist);
                    velocity[3*i+1] += dt*diffy/(dist*dist*dist);
                    velocity[3*i+2] += dt*diffz/(dist*dist*dist);
                }
            }
            double stop_compute = MPI_Wtime();
            t_compute += stop_compute-start_compute;

            // make sure transactions have finished
            MPI_Waitall(blocking?1:2, requests, statuses);
            // put recieved data into otherposition for next iteration
            std::swap(otherposition, recvbuf);
        }

        { // process last block of data and update positions
            double start_compute = MPI_Wtime();

            for(int i = 0;i<Nlocal;++i) {
                for(int j = 0;j<Nlocal;++j) {
                    float diffx = otherposition[3*i+0]-position[3*j+0];
                    float diffy = otherposition[3*i+1]-position[3*j+1];
                    float diffz = otherposition[3*i+2]-position[3*j+2];
                    float dist = std::sqrt(diffx*diffx + diffy*diffy + diffz*diffz)+0.001;
                    velocity[3*i+0] += dt*diffx/(dist*dist*dist);
                    velocity[3*i+1] += dt*diffy/(dist*dist*dist);
                    velocity[3*i+2] += dt*diffz/(dist*dist*dist);
                }
                position[3*i+0] += dt*velocity[3*i+0];
                position[3*i+1] += dt*velocity[3*i+1];
                position[3*i+2] += dt*velocity[3*i+2];
            }

            double stop_compute = MPI_Wtime();
            t_compute += stop_compute-start_compute;
        }
    }
    double stop = MPI_Wtime();
    double t_total = stop-start;

    // reduce times to node 0
    double t[2], avg[2];
    t[0] = t_total;
    t[1] = t_compute;
    MPI_Reduce(t, avg, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    avg[0] /= processes;
    avg[1] /= processes;

    if(rank == 0) {
        std::cout << "total:   " << avg[0]*1.e3 << "ms" << std::endl;
        std::cout << "compute: " << avg[1]*1.e3 << "ms" << std::endl;
    }
    // shutdown MPI
    MPI_Finalize();

    return 0;
}
