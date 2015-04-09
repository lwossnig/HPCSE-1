#include "timer.hpp"
#include <iostream>

float compute(unsigned long N, float mul,float sum,unsigned int nthreads) {
    float a=1.0,b=1.1,c=1.2,d=1.3,e=1.4,f=1.5,g=1.6,h=1.7,i=1.8,j=1.9; // assign different values to avoid compiler optimizations
    unsigned long loops = N/10;// loop contains 10 flops
#pragma omp parallel for num_threads(nthreads) reduction(*:a) reduction(+:b) reduction(*:c) reduction(+:d) reduction(*:e) reduction(+:f) reduction(*:g) reduction(+:h) reduction(*:i) reduction(+:j)
    for (unsigned long z = 0; z < loops; z++) {
        a*=mul; b+=sum; c*=mul; d+=sum; e*=mul; f+=sum; g*=mul; h+=sum; i*=mul; j+=sum;
    }
    return a+b+c+d+e+f+g+h+i+j;
}

 int main(int argc, char* argv[]){
    unsigned long n;
    unsigned int nthreads;
    switch(argc) {
    case 3:
        sscanf(argv[1],"%li",&n);
        sscanf(argv[2],"%i",&nthreads);
        break;
    case 1:                     // defaults
        n=1000;
        nthreads=1;
        break;
    default:
        std::cout << "Usage: " << argv[0] << " <number of operations (millions)> <num threads>" << std::endl;
        return 0;
        break;
    }
    n*=1000000;                 // millions of operations
    timer t;
    unsigned int N=10;
    float res,                 // return and use the value to avoid compiler optimizations
        avg=0;
    for (unsigned int i = 0; i < N; i++) { // compute time average
        t.start();
        res=compute(n,1.1,2.1,nthreads);     // multiply by a num !=1 to avoid compiler optimizations
        t.stop();
        avg+=t.get_timing();
    }
    avg=avg/N;
    std::cout << "Computed " << n/1000000 << " millions of operations in " <<  avg << " seconds: GFlops/s: " << n/avg/1e9<< " - res: " << res << std::endl;
}
