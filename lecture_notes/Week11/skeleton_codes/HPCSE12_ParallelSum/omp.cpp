// Skeleton code for HPCSE Exam, 18.12.2012
// Profs. P. Koumoutsakos and M. Troyer
// Question 5a)

#include <vector>
#include <numeric>
#include <iostream>
#include <omp.h>

int main( int argc, char** argv )
{
    // vector size
    const int N = 1600000;

    // initialize vectors
    std::vector<float> x(N,-1.2), y(N,3.4), z(N);
    
    
    // DO THE SUM z = x + y
#pragma omp parallel for
    for(int i=0; i<N; i++)
	    z[i] = x[i] + y[i];

    
    // print result checksum
    std::cout << std::accumulate(z.begin(), z.end(), 0.) << std::endl;
}
