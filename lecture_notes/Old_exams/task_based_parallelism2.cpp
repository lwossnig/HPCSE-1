#include <iostream>
#include <thread>
#include <chrono>
#include <unistd.h>
#include <vector>
#include <mutex>

int alpha () {	sleep(1);	return 10;  }
int beta  () {	sleep(1);	return 20;  }
int gamma () {	sleep(1);	return 50;  }
int delta () {	sleep(1);	return 80;  }
int epsilon ( int x, int y, int z )  { sleep(1);     return x*y+z;  }
int zeta ( int x, int y, int z )     { sleep(1);     return x-y+z;  }
int ita ( int x, int y )	     { sleep(1);     return x/y;    }


void wait(std::pair<int, std::mutex>& m, int i){
	std::unique_lock<std::mutex> l(m.second);
	l.unlock();
	while(true){
		l.lock();
		if(m.first >= i)
			break;
		l.unlock();
	}
}

int main()
{
	int A , B, C, D, E, F, G;
	std::vector<std::thread> threads(7);

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
	std::pair<int, std::mutex> count;
	count.first=0;
	std::unique_lock<std::mutex> l(count.second);
	l.unlock();
	t0 = std::chrono::high_resolution_clock::now();
	std::cout << "First threads running!" << std::endl;

	threads[0] = std::thread(
			[&](){ 
			A = alpha ();
			l.lock();
			count.first++;
			l.unlock();
			});
	
	threads[1] = std::thread(
			[&](){ 
			B = beta  ();
			l.lock();
			count.first++;
			l.unlock();
			});

	threads[2] = std::thread(
			[&](){ 
			C = gamma (); 
			l.lock();
			count.first++;
			l.unlock();
			});

	threads[3] = std::thread(
			[&](){ 
			D = delta ();
			l.lock();
			count.first++;
			l.unlock();
			});
	

	
	threads[4] = std::thread(
			[&](){ 
			wait(count, 4);
			E = epsilon (A, B, C);
			l.lock();
			count.first++;
			l.unlock();
			});


	threads[5] = std::thread(
			[&](){ 
			wait(count, 4);
			F = zeta    (B, C, D); 
			l.lock();
			count.first++;
			l.unlock();
			});
	

	threads[6] = std::thread(
			[&](){
			wait(count, 6);
			G = ita	(E, F); 
			});

	for(size_t i =0; i<7; ++i)
		threads[i].join();


	t1 = std::chrono::high_resolution_clock::now();
	elapsed_time = std::chrono::duration<double>(t1-t0).count() ;

	std::cout << "G threaded = " << G << std::endl;
	std::cout << "Elapsed threaded = " << elapsed_time << std::endl;
	
	return 0;
}

