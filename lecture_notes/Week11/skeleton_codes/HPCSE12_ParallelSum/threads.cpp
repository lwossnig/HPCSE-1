// Skeleton code for HPCSE Exam, 18.12.2012
// Profs. P. Koumoutsakos and M. Troyer
// Question 5b)

#include <vector>
#include <numeric>
#include <iostream>
#include <thread>

int main( int argc, char** argv )
{
    // vector size
    const int N = 1600000;

    // initialize vectors
    std::vector<float> x(N,-1.2), y(N,3.4), z(N);
    std::vector<std::thread> threads(4);
    int chunksize = N/4;
    
    // DO THE SUM z = x + y using 4 threads
    for(size_t i=0; i<4; ++i){
	    threads[i] = std::thread(
			    [&,i]() {
			    for(size_t j = i*chunksize; j < (i+1)*chunksize; ++j)
				    z[j] = x[j] + y[j];
			    });
    }
    for(size_t i =0;i <4; ++i)
	    threads[i].join();
    
    // print result checksum
    std::cout << std::accumulate(z.begin(), z.end(), 0.) << std::endl;
}

