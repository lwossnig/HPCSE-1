#include <iostream>
#include <omp.h>
#include <chrono>
#include <unistd.h>
#include <vector>

int alpha () {	sleep(1);	return 10;  }
int beta  () {	sleep(1);	return 20;  }
int gamma () {	sleep(1);	return 50;  }
int delta () {	sleep(1);	return 80;  }
int epsilon ( int x, int y, int z )  { sleep(1);     return x*y+z;  }
int zeta ( int x, int y, int z )     { sleep(1);     return x-y+z;  }
int ita ( int x, int y )	     { sleep(1);     return x/y;    }


int main()
{
	int A , B, C, D, E, F, G;

	std::chrono::time_point<std::chrono::high_resolution_clock> t0, t1;
	t0 = std::chrono::high_resolution_clock::now();

	A = alpha();
	B = beta ();
	C = gamma();
	D = delta();
	E = epsilon(A, B, C);
	F = zeta   (B, C, D);
	G = ita	   (E, F);


	t1 = std::chrono::high_resolution_clock::now();
	double elapsed_time = std::chrono::duration<double>(t1-t0).count() ;

	std::cout << "G = " << G << std::endl;
	std::cout << "Elapsed = " << elapsed_time << std::endl;

	/***** HERE THREADED *****/
	A=0;B=0;C=0;D=0;E=0;F=0;G=0;
	t0 = std::chrono::high_resolution_clock::now();

#pragma omp parallel num_threads(4)
	{
		if(omp_get_thread_num()==0)
			A = alpha ();
		if(omp_get_thread_num()==1)
			B = beta  ();
		if(omp_get_thread_num()==2)
			C = gamma (); 
		if(omp_get_thread_num()==3)
			D = delta ();

#pragma omp barrier
		if(omp_get_thread_num()==0)
		{
			E = epsilon (A, B, C);
		}
		if(omp_get_thread_num()==1)
		{
			F = zeta    (B, C, D); 
		}
#pragma omp barrier
#pragma omp single
		{
			G = ita	(E, F); 
		}
	}
	t1 = std::chrono::high_resolution_clock::now();
	elapsed_time = std::chrono::duration<double>(t1-t0).count() ;

	std::cout << "G threaded = " << G << std::endl;
	std::cout << "Elapsed threaded = " << elapsed_time << std::endl;
	
	return 0;
}

