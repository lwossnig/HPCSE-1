/* The following code was derived and developed following "Parallel Numerical Solution of
 * 2-D Heat Equation", Verena Horak, Peter Gruber - Department of Scientific Computing,
 * Univ. of Salzburg - Parallel Numerics '05, 47-56 Chapter 3: Differential Equations 
 * by M.Vaijtersic, R Trobec, P. Zinterhof, Aa. Uhl,
 * ISBN 961-6303-67-8
 * The fundamental ideas of this algorithm are derived in there, yet the code itself was 
 * written and applied to the exercise sheet by Leonard Wossnig and the indidual derivatives
 * were again calculated specifically and by L.W. himself again.
 */

#include<iostream>
#include<cmath>
#include<vector>
#include<iomanip>
#include<cassert>
using namespace std;

int main(int argc, char *argv[])
{
	// Building a matrix where the rows are the correponding y values and the colum are the x values where the ones	  
	double D, time, delta_t, delta_s;
	unsigned int nthreads = 1;
	int steps;
	int time_steps = 1000;
	time = 1.d;
	D = 1.d;
	steps = 12; 
	if (argc > 1) time_steps = atoi(argv[1]);
	if (argc > 2) steps = atoi(argv[2]);
	if (argc > 3) nthreads = atoi(argv[3]);
	if (argc > 4) time = atoi(argv[4]);
	assert(steps%4==0);
	delta_t = 1./time_steps;
	delta_s = 2./(steps-1); // since interval [-1,1]
      	unsigned int N = (steps)*(steps);
	double alpha = delta_t * D/(delta_s*delta_s); // factor alpha for iteration

	// initialize vectors with steps * steps entries = grid for calculation roh(i,j)
	vector<double> u_old(N,0);
	vector<double> u_new(N,0);
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
    for(int k=0;k<time;k++)
      {
	// iteration through the lattice and applying the algorithm
	for(int i=1;i<(steps-1);i++)
	{
		for(int j=1;j<(steps-1);j++)
		{
			u_new[ i * steps + j] = u_new[ i * steps + j] + alpha*( u_old[i * steps + j + 1] - 4 * u_old[ i * steps + j] + u_old[i * steps + j - 1] + u_old[i * steps + j + steps] + u_old[i * steps + j - steps]);
		}
	}


	u_old.swap(u_new);
      }  	
      
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


/* In case for plotting the matrix to the terminal
   	for(int i=0;i<steps;i++)
	{
		for(int j =0;j<steps;j++)
		{
			cout << u_new[i*steps+j] << " ";
		}
		cout << endl;
	}
*/	
