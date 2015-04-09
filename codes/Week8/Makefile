
run: vortex_mpi
	mpirun.mpich -n 8 ./vortex_mpi

vortex_mpi: Vortex_mpi.cpp
	mpic++.mpich -std=c++11 -o vortex_mpi Vortex_mpi.cpp

