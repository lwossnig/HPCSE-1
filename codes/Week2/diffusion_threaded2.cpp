/* The following code was derived and developed following "Parallel Numerical Solution of
 * 2-D Heat Equation", Verena Horak, Peter Gruber - Department of Scientific Computing,
 * Univ. of Salzburg - Parallel Numerics '05, 47-56 Chapter 3: Differential Equations 
 * by M.Vaijtersic, R Trobec, P. Zinterhof, Aa. Uhl,
 * ISBN 961-6303-67-8
 * The fundamental ideas of this algorithm are derived in there, yet the code itself was 
 * written and applied to the exercise sheet by Leonard Wossnig and the indidual derivatives
 * were again calculated specifically and by L.W. himself again.
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <cassert>
#include <thread>
#include "timer.hpp"
#include <mutex>
#include "syncedio.cpp"

using namespace std;

int main(int argc, char *argv[])
{
	// Building a matrix where the rows are the correponding y values and the colum are the x values where the ones	  
	double D, time, delta_t, delta_s;
	unsigned int nthreads = 1; //number of threads
	int steps; //number of steps
	int time_steps = 10000;
	time = 1.d; //number of iterations
	D = 1.d; // Diffusion constant
	steps = 12; //# of steps 
	//Input while running?
	if (argc > 1) time_steps = atoi(argv[1]);
	if (argc > 2) steps = atoi(argv[2]);
	if (argc > 3) nthreads = atoi(argv[3]);
	if (argc > 4) time = atoi(argv[4]);
	// step size delta_x=delta_y = (1-(-1))/steps = 2/steps
	delta_t = 1.d/time_steps;
	delta_s = 2.d/(steps-1); // since interval [-1,1]
	int count_lock(nthreads);
	//cout << "steps: " << steps <<" - delta_t: " << delta_t << " -  nthreads: " << nthreads << " - time: " << time <<  endl;


	// ASSERT if grid is suitable (means if gitsize/4 is even), since we implemented the grid so that we just use the ratio [-1,1]
	// and the corresponding initial values |x,y|<0,5 = 1 are on grid just the inner values of the grip position: 1/4 * steps - 3/4 * steps
	assert(steps%4==0);


	unsigned int N = (steps)*(steps); // Arrazsize to work on
	assert(steps%nthreads==0); // check if Number of iterations can be divided equally through the n# of threads
	unsigned int n = (steps-2)/nthreads;
	double alpha = delta_t * D/(delta_s*delta_s); // coefficient from the iteration rule

	// initialize vectors with steps * steps entries = grid for calculation roh(i,j)
	vector<double> u_old(N,0);
	vector<double> u_new(N,0);
	//vector<bool> lock(nthreads,false); //initialize the "control" vector with all stati: false
	//vector<bool> status(nthreads,true);
	
	vector<thread> threads(nthreads);
	
	vector<pair<bool, mutex>> status(nthreads);
	for(auto& it : status) it.first =true;
	vector<pair<bool,mutex>> lock(nthreads);
	for(auto& it2 : lock) it2.first = false;
	//mutex status_mutex;
	//mutex lock_mutex;
	mutex array_lock;


	// Null initialisation so not neccessary to set the boundary conditions to zero
	
	// initial values: u(x,y,t) = 1 f.a. |x,y|<1/2
	for(int i=(0.25*steps);i<(0.75*steps);i++)
	{
		for(int j=(0.25*steps);j<(0.75*steps);j++)
		{
			u_old[i*steps+j] = 1.;
			u_new[i*steps+j] = 1.;
		}
	}

	// By dropping the summation inside the iteration for i,j=0,n we fullfill the dirichlet boundary conditions	

  
	// timer for the calculation
	timer t;
	t.start();

	
	
	// spawn threads
	for(unsigned int thr = 0; thr < nthreads; thr++)
	{
		
		threads[thr] = thread([&,thr]()
				{
				int wait = 1;
				for(int t=0;t<time;t++)
				{
				//sync(cout) << "thread "<< thr <<" running " << t <<" th time!\n";
				wait = 1;
				while(wait==1){
					      status[thr].second.lock();
					      if(status[thr].first==true) // if main didn't set status of thr to true it waits
							 {wait=0;}  // else it may go on (change wait to zero and escape while loop
					      status[thr].second.unlock();
					      }
				status[thr].second.lock();
				status[thr].first = false;
				status[thr].second.unlock();
				//sync(cout) << "status changed to " << status[thr].first << endl;
				
				
				//call lambda function
				// iteration through the lattice and applying the algorithm
				for(int i=1+thr*n;i<(thr+1)*n;i++) // i<(thr+1)*n
				    {
					//assure the boundaries are not touched!
					for(int j=1;j<(steps-1);j++)
					  {
						u_new[i*steps+j]=u_new[i*steps+j] + alpha*(u_old[i*steps+j+1] - 4*u_old[i*steps+j] + u_old[i*steps+j-1]	+ u_old[i*steps+j+steps] + u_old[i*steps+j-steps]);
					  }
				     } 
				

				     //change manual lock which prevent main to work on
				     //sync(cout) <<"thread "<< thr << " progressing!\n";
				     lock[thr].second.lock();
				     lock[thr].first =true;
				     lock[thr].second.unlock();
				}
				});
	}
	
	// use barrier(nthreads) to wait till everybody updated, then swap!
	//prevent main to join threads before the added results
	for(int t=0;t<time;t++)
	{
	int help_var;
	while(count_lock>0)
	{
	    help_var=count_lock;
	    for(int k=0;k<nthreads;k++)
	 	  {
			  lock[k].second.lock();
	 	          if(lock[k].first == true){
				  count_lock--; //if the thread did it's job and set lock to true, lock the thread to go on (status = false)
			  }
			  lock[k].second.unlock();
	 	   }
	    if(count_lock!=0)count_lock=help_var;
	 }
	//sync(cout) << "Main progressing in " << t<<" th round!\n";
	
	// make new time step
	//sync(cout)<<"swapping\n";
	array_lock.lock();
	u_old.swap(u_new);
	array_lock.unlock();
	//sync(cout) << "swapped!\n";
	// set again all to false to begin new loop
	count_lock = nthreads;
		for(int i=0;i<nthreads;i++)
	 	  {
			  lock[i].second.lock();
		          lock[i].first=false;
			  lock[i].second.unlock();
	 	  }
	
	// free threads to go on
	for(int t=0;t<nthreads;t++){
		status[t].second.lock();
		status[t].first=true;
		status[t].second.unlock();
	}
		//sync(cout) << "All stati set to true!\n";
	}

	// End threads
	for(thread& thr: threads)
		  thr.join();
	
	t.stop();
	
	
	//cout << "Time for " << nthreads << " threads: " << t.get_timing() << "seconds" << endl;

	// Output resulting matrix (points) 
    	for(int i=0;i<steps;i++)
	{
		for(int j =0;j<steps;j++)
		{
			cout << setw(7) << setprecision(3) << u_new[i*steps+j] << "  ";
		}
		cout << endl;
	}
 

  return 0;
}

