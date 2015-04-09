#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cassert>

using namespace std;

class PSE2D
{
private:
    const int support[2] = {-10,11};
    const double eps;
    
public:
    PSE2D(int N) : eps(2./N) {}
    
    inline const int getStart() const
    {
        return support[0];
    }
    
    inline const int getEnd() const
    {
        return support[1];
    }
    
    inline const double w(double dx, double dy) const
    {
        //PSE kernel evaluation
	const double dist = dx * dx + dy * dy;
	const double d2	  = dist/(eps*eps);
	return 16./(M_PI*M_PI) /( d2*d2*d2*d2 + 1. );
    }
};


template<class Kernel>
class DiffusionPSE
{
private:
    // data structures for position (x) and concentrations (c)
    // these array contain ghosts in addition to the domain of each rank
    vector<double> cPing, cPong;
    
    // PSE kernel
    Kernel k;
    
    // Grid spacing, problem size
    double dh;
    const int N;
    
    // MPI stuff
    int rank, size, Mx, My;
    int rank_x, rank_y;
    int start_x, start_y, end_x, end_y, worksize_x, worksize_y, workspace_x, workspace_y;
    int bufsizeCorner, bufsizeFace_x, bufsizeFace_y;

    // MPI communication buffers
    // Buffers for all 8 neighbours (note, here only pointers for later dynamic buffers!)
    double * bufSendN;
    double * bufSendS;
    double * bufSendW;
    double * bufSendE;
    double * bufSendNW;
    double * bufSendNE;
    double * bufSendSW;
    double * bufSendSE;
    double * bufRecvN;
    double * bufRecvS;
    double * bufRecvW;
    double * bufRecvE;
    double * bufRecvNW;
    double * bufRecvNE;
    double * bufRecvSW;
    double * bufRecvSE;

    double analytical(double x, double y, double t)
    {
        return sin(2.*M_PI*x)*sin(2.*M_PI*y)*exp(-8*M_PI*t);
    }
    
public:
    DiffusionPSE(   int rank, 
		    int size, 
		    int N, 
		    int Mx, 
		    int My) 
	    : rank(rank), size(size), dh(1./(double)N), k(N), N(N), Mx(Mx), My(My)
    {
	    if (N % Mx != 0)
	    {
		    if (rank==0) cout << "Problem size is not divisible by number of ranks in x-direction!\nAborting now...\n";
		    abort();
	    }
	    if (N % My != 0)
	    {
		    if (rank==0) cout << "Problem size is not divisible by number of ranks in y-direction!\nAborting now...\n";
		    abort();
	    }

	    /// initialize MPI
	    
	    rank_x = rank % Mx;
	    rank_y = rank / Mx;
	    
	    worksize_x = N / Mx;
	    worksize_y = N / My;
	    start_x = -k.getStart();
	    start_y = -k.getStart();
	    end_x = worksize_x - k.getStart();
	    end_y = worksize_y - k.getStart();
	    workspace_x = worksize_x - k.getStart() + k.getEnd() -1;
	    workspace_y = worksize_y - k.getStart() + k.getEnd() -1;

	    cPing = vector<double>(workspace_x * workspace_y,0);
	    cPong = vector<double>(workspace_x * workspace_y,0);

	    assert(-k.getStart() == k.getEnd()-1);
	    bufsizeCorner = start_x * start_y;
	    bufsizeFace_x = start_y * worksize_x;
	    bufsizeFace_y = start_x * worksize_y;


	    bufSendN = new double[bufsizeFace_x];
	    bufSendS = new double[bufsizeFace_x];
	    bufSendW = new double[bufsizeFace_y];
	    bufSendE = new double[bufsizeFace_y];
	    bufSendNW = new double[bufsizeCorner];
	    bufSendNE = new double[bufsizeCorner];
	    bufSendSW = new double[bufsizeCorner];
	    bufSendSE = new double[bufsizeCorner];
	    bufRecvN = new double[bufsizeFace_x];
	    bufRecvS = new double[bufsizeFace_x];
	    bufRecvW = new double[bufsizeFace_y];
	    bufRecvE = new double[bufsizeFace_y];
	    bufRecvNW = new double[bufsizeCorner];
	    bufRecvNE = new double[bufsizeCorner];
	    bufRecvSW = new double[bufsizeCorner];
	    bufRecvSE = new double[bufsizeCorner];

	    // D=1, dh=1/N
	    // eps=2*dx=2*dy=2*dh
	    // dt= dh*dh/(2*D)
	    coeff = 0.5 / double(N*N);

	    // additional factor D / eps^2 * vol = 1 / (2dh)^2 * dh^2
	    coeff *= 0.25;
    }
    
    ~DiffusionPSE()
    {
	    delete []	bufSendN  ;
	    delete []	bufSendS  ;
	    delete []	bufSendW  ;
	    delete []	bufSendE  ;
	    delete []	bufSendNW ;
	    delete []	bufSendNE ;
	    delete []	bufSendSW ;
	    delete []	bufSendSE ;
	    delete []	bufRecvN  ;
	    delete []	bufRecvS  ;
	    delete []	bufRecvW  ;
	    delete []	bufRecvE  ;
	    delete []	bufRecvNW ;
	    delete []	bufRecvNE ;
	    delete []	bufRecvSW ;
	    delete []	bufRecvSE ;
    }
    
    void setIC()
    {
	    for	(int iy=start_y; iy<end_y; iy++)
		    for	(int ix=start_x; ix<end_x; ix++)
		    {
			    const int gix = rank_x * worksize_x + ix + k.getStart();
			    const int giy = rank_y * worksize_y + iy + k.getStart();
			    cPing [ix + iy * workspace_x] = analytical	( ( double)gix/(double)N, (double)giy/(double)N, 0  );
		    }
	    double localC =0, totalC;
	    for(int iy=start_y; iy<end_y; iy++)
		    for(int ix=start_x; ix<end_x; ix++)
			    localC += cPing[ix + iy * workspace_x];

	    MPI_Reduce(&localC, &totalC, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	    if (rank==0)
		    cout << "\nTotal concentration at IC is " << totalC << endl;
    }

    void simulate()
    {
        // time loop
        const int nsteps = 100;
        const int supportSize = -k.getStart();
        for (int t=0; t<nsteps; t++)
        {
           // pack corners
	   for	(int iy=0; iy < start_y; iy++)
		   for	(int ix=0; ix < start_x; ix++)
		   {
			   const int xs = start_x + ix;
			   const int ys = start_y + iy;
			   const int xe = end_x	- start_x + ix;
			   const int ye = end_y - start_y + iy;
			   bufSendNE [ix + iy * start_x] = cPing [xe + workspace_x * ye];
			   bufSendNW [ix + iy * start_x] = cPing [xs + workspace_x * ye];
			   bufSendSE [ix + iy * start_x] = cPing [xe + workspace_x * ye];
			   bufSendSW [ix + iy * start_x] = cPing [xs + workspace_x * ye];
		   }




            // swap the two data structures s.t. the current solution is always in cPing
            cPing.swap(cPong);
        }
    }
};

int main(int argc, char *argv[])
{
    //=========================
    //      Initialize MPI
    //=========================
    MPI_Init(&argc, &argv);
    
    // set MPI configuration
    Mx = 2;
    My = 2;
  
    //=========================
    //      get rank information
    //      (processID and #processes)
    //=========================
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != Mx*My)
    {
        if (rank==0) cout << "Inconsistent number of MPI Processes\n";
        abort();
    }
    printf("Hello from process %d of %d\n", rank, size);
    
    // Problem size
    int N = 120;
    
    // setup simulation
    DiffusionPSE<PSE2D> diffusion(rank,size,N,Mx,My);
    diffusion.setIC();
    
    // run simulation
    diffusion.simulate();
    
    
    //=========================
    //      cleanup MPI environment
    //=========================
    MPI_Finalize();
    
    return 0;
}
