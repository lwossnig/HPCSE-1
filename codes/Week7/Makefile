run: diffusion2d_mpi
	mpirun.mpich -n 4 ./diffusion2d_mpi 1 1 128 0.00001

diffusion2d_mpi:
	mpic++.mpich -std=c++11 -o diffusion2d_mpi diffusion2d_mpi.cpp


