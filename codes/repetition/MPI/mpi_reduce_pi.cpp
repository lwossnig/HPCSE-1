#include <iostream>
#include <cmath>
#include <iomanip>
#include <mpi.h>

int main  (int argc, char* argv[])  {
      MPI_Init(&argc, &argv);
      int rank, size;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      long double sum =0.;
      long double localsum =0.;
      unsigned long const nterms = 100000000;
      long double const step = (nterms+0.5l) / size;

      unsigned long start = rank * step;
      unsigned long stop = (rank+1) * step;
      for (std::size_t t = start; t < stop; ++t)
	      localsum += (1.0 -2* (t % 2)) / (2*t + 1);

      MPI_Reduce(&localsum, &sum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      if  (rank==0)
	      std::cout << "pi = " << std::setprecision(18) << 4. * sum << std::endl;

      MPI_Finalize();
      return 0;
}

