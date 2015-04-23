1) Load modules
module load gcc
module load open_mpi

2) Compile
mpicc -o example1 example1.c 

3) Submit script
bsub -n 48 < script


