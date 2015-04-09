//
//  main.cpp
//
//  2D Diffusion with PSE using MPI (ghosts version)
//
//  Created by Christian Conti on 12/13/13.
//
//

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cassert>

using namespace std;

// mod function for integers
inline int mod (const int &a, const int &b) {
    int m = a % b;
    return (m >= 0) ? m : m + b;
}

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
        const double dist = dx*dx + dy*dy;
        const double d2 = dist/(eps*eps);
        return 16./(M_PI*M_PI)/(d2*d2*d2*d2 + 1.);
    }
};

template<class Kernel>
class DiffusionPSE
{
private:
    // data structures for position (x) and concentrations (c)
    // these array contain ghosts in addition to the domain of each rank
    vector<double> cPing, cPong;
    
    // Diffusion kernel
    Kernel k;
    
    // Grid spacing, problem size
    double coeff, dh;
    const int N;
    
    // MPI stuff
    int rank, size, Mx, My;
    int rank_x, rank_y;
    int start_x, start_y, end_x, end_y, worksize_x, worksize_y, workspace_x, workspace_y;
    int bufsizeCorner, bufsizeFace_x, bufsizeFace_y;
    
    // MPI communication buffers
    double * bufSendN;
    double * bufSendS;
    double * bufSendW;
    double * bufSendE;
    double * bufSendNE;
    double * bufSendNW;
    double * bufSendSE;
    double * bufSendSW;
    double * bufRecvN;
    double * bufRecvS;
    double * bufRecvW;
    double * bufRecvE;
    double * bufRecvNE;
    double * bufRecvNW;
    double * bufRecvSE;
    double * bufRecvSW;
    
    double analytical(double x, double y, double t)
    {
        return sin(2.*M_PI*x)*sin(2.*M_PI*y)*exp(-8*M_PI*t);
    }
    
public:
    DiffusionPSE(int rank, int size, int N, int Mx, int My) : rank(rank), size(size), dh(1./(double)N), k(N), N(N), Mx(Mx), My(My)
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
      
        rank_x = rank % Mx;
        rank_y = rank / Mx;
        
        worksize_x = N / Mx;
        worksize_y = N / My;
        start_x = -k.getStart();
        start_y = -k.getStart();
        end_x = worksize_x - k.getStart();
        end_y = worksize_y - k.getStart();
        workspace_x = worksize_x - k.getStart() + k.getEnd()-1;
        workspace_y = worksize_y - k.getStart() + k.getEnd()-1;
        
        cPing = vector<double>(workspace_x * workspace_y,0);
        cPong = vector<double>(workspace_x * workspace_y,0);
        
        assert(-k.getStart() == k.getEnd()-1);
        
        bufsizeCorner = start_x * start_y;
        bufsizeFace_x = start_y * worksize_x;
        bufsizeFace_y = start_x * worksize_y;
        
        bufSendN  = new double[bufsizeFace_x];
        bufSendS  = new double[bufsizeFace_x];
        bufSendW  = new double[bufsizeFace_y];
        bufSendE  = new double[bufsizeFace_y];
        bufSendNE = new double[bufsizeCorner];
        bufSendNW = new double[bufsizeCorner];
        bufSendSE = new double[bufsizeCorner];
        bufSendSW = new double[bufsizeCorner];
        bufRecvN  = new double[bufsizeFace_x];
        bufRecvS  = new double[bufsizeFace_x];
        bufRecvW  = new double[bufsizeFace_y];
        bufRecvE  = new double[bufsizeFace_y];
        bufRecvNE = new double[bufsizeCorner];
        bufRecvNW = new double[bufsizeCorner];
        bufRecvSE = new double[bufsizeCorner];
        bufRecvSW = new double[bufsizeCorner];
        
        // D=1, 1/dh=N
        // dt = D / 2 * dh^2
        coeff = 0.5 / double(N*N);
      
        // additional factor D / eps^2 * vol = 1 / (2dh)^2 * dh^2
        coeff *= 0.25;
    }
    
    ~DiffusionPSE()
    {
        delete [] bufSendN;
        delete [] bufSendS;
        delete [] bufSendW;
        delete [] bufSendE;
        delete [] bufSendNW;
        delete [] bufSendNE;
        delete [] bufSendSW;
        delete [] bufSendSE;
        delete [] bufRecvN;
        delete [] bufRecvS;
        delete [] bufRecvW;
        delete [] bufRecvE;
        delete [] bufRecvNW;
        delete [] bufRecvNE;
        delete [] bufRecvSW;
        delete [] bufRecvSE;
    }
    
    void setIC()
    {
        for (int iy=start_y; iy<end_y; iy++)
            for (int ix=start_x; ix<end_x; ix++)
            {
                const int gix = rank_x * worksize_x + ix + k.getStart();
                const int giy = rank_y * worksize_y + iy + k.getStart();
                cPing [ix + iy * workspace_x] = analytical ( (double)gix/(double)N, (double)giy/(double)N, 0 );
            }
        
        double localC = 0, totalC;
        for (int iy=start_y; iy<end_y; iy++)
            for (int ix=start_x; ix<end_x; ix++)
                localC += cPing [ix + iy * workspace_x];
        
        MPI_Reduce(&localC, &totalC, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if (rank==0)
            cout << "\nTotal concentration at IC is " << totalC << endl;
    }
    
    void simulate()
    {
        // time loop
        const int nsteps = 100;
        for (int t=0; t<nsteps; t++)
        {
            // pack corners
            for (int iy=0; iy<start_y; iy++)
                for (int ix=0; ix<start_x; ix++)
                {
                    const int xs = start_x + ix;
                    const int ys = start_y + iy;
                    const int xe = end_x - start_x + ix;
                    const int ye = end_y - start_y + iy;
                    bufSendNE [ix + iy * start_x] = cPing [xe + workspace_x * ye];
                    bufSendNW [ix + iy * start_x] = cPing [xs + workspace_x * ye];
                    bufSendSE [ix + iy * start_x] = cPing [xe + workspace_x * ys];
                    bufSendSW [ix + iy * start_x] = cPing [xs + workspace_x * ys];
                }
            
            // pack N and S faces
            for (int iy=0; iy<start_y; iy++)
                for (int ix=0; ix<worksize_x; ix++)
                {
                    const int x  = start_x + ix;
                    const int yS = start_y + iy;
                    const int yN = end_y - start_y + iy;
                    bufSendN [ix + worksize_x * iy] = cPing [x + workspace_x * yN];
                    bufSendS [ix + worksize_x * iy] = cPing [x + workspace_x * yS];
                }
            
            // pack E and W faces
            for (int iy=0; iy<worksize_y; iy++)
                for (int ix=0; ix<start_x; ix++)
                {
                    const int y  = start_y + iy;
                    const int xW = start_x + ix;
                    const int xE = end_x - start_x + ix;
                    bufSendW [ix + start_x * iy] = cPing [xW + workspace_x * y];
                    bufSendE [ix + start_x * iy] = cPing [xE + workspace_x * y];
                }
            
            MPI_Request reqs[16];
            
            const int TAG_S = 0;
            const int TAG_N = 1;
            const int TAG_W = 2;
            const int TAG_E = 3;
            
            const int TAG_SW = 4;
            const int TAG_SE = 5;
            const int TAG_NW = 6;
            const int TAG_NE = 7;
            
            const int rank_S = rank_x + Mx * mod (rank_y-1, My);
            const int rank_N = rank_x + Mx * mod (rank_y+1, My);
            const int rank_W = mod (rank_x-1, Mx) + Mx * rank_y;
            const int rank_E = mod (rank_x+1, Mx) + Mx * rank_y;
            
            const int rank_SW = mod (rank_x-1, Mx) + Mx * mod (rank_y-1, My);
            const int rank_SE = mod (rank_x+1, Mx) + Mx * mod (rank_y-1, My);
            const int rank_NW = mod (rank_x-1, Mx) + Mx * mod (rank_y+1, My);
            const int rank_NE = mod (rank_x+1, Mx) + Mx * mod (rank_y+1, My);
            
            MPI_Irecv (bufRecvS,  bufsizeFace_x, MPI_DOUBLE, rank_S,  TAG_S,  MPI_COMM_WORLD, &reqs[0] );
            MPI_Irecv (bufRecvN,  bufsizeFace_x, MPI_DOUBLE, rank_N,  TAG_N,  MPI_COMM_WORLD, &reqs[1] );
            MPI_Irecv (bufRecvW,  bufsizeFace_y, MPI_DOUBLE, rank_W,  TAG_W,  MPI_COMM_WORLD, &reqs[2] );
            MPI_Irecv (bufRecvE,  bufsizeFace_y, MPI_DOUBLE, rank_E,  TAG_E,  MPI_COMM_WORLD, &reqs[3] );
            
            MPI_Irecv (bufRecvSW, bufsizeCorner, MPI_DOUBLE, rank_SW, TAG_SW, MPI_COMM_WORLD, &reqs[4] );
            MPI_Irecv (bufRecvSE, bufsizeCorner, MPI_DOUBLE, rank_SE, TAG_SE, MPI_COMM_WORLD, &reqs[5] );
            MPI_Irecv (bufRecvNW, bufsizeCorner, MPI_DOUBLE, rank_NW, TAG_NW, MPI_COMM_WORLD, &reqs[6] );
            MPI_Irecv (bufRecvNE, bufsizeCorner, MPI_DOUBLE, rank_NE, TAG_NE, MPI_COMM_WORLD, &reqs[7] );
            
            MPI_Isend (bufSendS,  bufsizeFace_x, MPI_DOUBLE, rank_S,  TAG_N,  MPI_COMM_WORLD, &reqs[8] );
            MPI_Isend (bufSendN,  bufsizeFace_x, MPI_DOUBLE, rank_N,  TAG_S,  MPI_COMM_WORLD, &reqs[9] );
            MPI_Isend (bufSendW,  bufsizeFace_y, MPI_DOUBLE, rank_W,  TAG_E,  MPI_COMM_WORLD, &reqs[10]);
            MPI_Isend (bufSendE,  bufsizeFace_y, MPI_DOUBLE, rank_E,  TAG_W,  MPI_COMM_WORLD, &reqs[11]);
            
            MPI_Isend (bufSendSW, bufsizeCorner, MPI_DOUBLE, rank_SW, TAG_NE, MPI_COMM_WORLD, &reqs[12]);
            MPI_Isend (bufSendSE, bufsizeCorner, MPI_DOUBLE, rank_SE, TAG_NW, MPI_COMM_WORLD, &reqs[13]);
            MPI_Isend (bufSendNW, bufsizeCorner, MPI_DOUBLE, rank_NW, TAG_SE, MPI_COMM_WORLD, &reqs[14]);
            MPI_Isend (bufSendNE, bufsizeCorner, MPI_DOUBLE, rank_NE, TAG_SW, MPI_COMM_WORLD, &reqs[15]);
            
            MPI_Waitall (16, reqs, MPI_STATUS_IGNORE);
            
            // unpack corners
            for (int iy=0; iy<start_y; iy++)
                for (int ix=0; ix<start_x; ix++)
                {
                    cPing [ix           + workspace_x * (end_y + iy)] = bufRecvNW [ix + iy * start_x];
                    cPing [(end_x + ix) + workspace_x * (end_y + iy)] = bufRecvNE [ix + iy * start_x];
                    cPing [ix           + workspace_x * iy]           = bufRecvSW [ix + iy * start_x];
                    cPing [(end_x + ix) + workspace_x * iy]           = bufRecvSE [ix + iy * start_x];
                }
            
            // unpack N and S faces
            for (int iy=0; iy<start_y; iy++)
                for (int ix=0; ix<worksize_x; ix++)
                {
                    const int x  = start_x + ix;
                    const int yS = iy;
                    const int yN = end_y + iy;
                    cPing [x + workspace_x * yN] = bufRecvN [ix + worksize_x * iy];
                    cPing [x + workspace_x * yS] = bufRecvS [ix + worksize_x * iy];
                }
            
            // unpack E and W faces
            for (int iy=0; iy<worksize_y; iy++)
                for (int ix=0; ix<start_x; ix++)
                {
                    const int y  = start_y + iy;
                    const int xW = ix;
                    const int xE = end_x + ix;
                    cPing [xW + workspace_x * y] = bufRecvW [ix + start_x * iy];
                    cPing [xE + workspace_x * y] = bufRecvE [ix + start_x * iy];
                }
            
            // apply PSE kernel
            for (int iy=start_y; iy<end_y; iy++)
                for (int ix=start_x; ix<end_x; ix++)
                {
                    // gather operation
                    double tmp = 0;
                    for (int sy=k.getStart(); sy<k.getEnd(); sy++)
                        for (int sx=k.getStart(); sx<k.getEnd(); sx++)
                            tmp += (cPing [ix+sx + workspace_x * (iy+sy)] - cPing [ix + workspace_x * iy]) * k.w((double)sx*dh, (double)sy*dh);
                    
                    cPong [ix + workspace_x * iy] = cPing [ix + workspace_x * iy] + tmp * coeff;
                }
            
            // swap the two data structures s.t. the current solution is always in cPing
            cPing.swap (cPong);
        }
        
        // check total concentration
        double localC = 0, totalC;
        for (int iy=start_y; iy<end_y; iy++)
            for (int ix=start_x; ix<end_x; ix++)
                localC += cPing [ix + iy * workspace_x];
        
        //=========================
        //      MPI_Reduce
        //=========================
        MPI_Reduce(&localC, &totalC, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if (rank==0)
            cout << "\nTotal concentration at end of simulation is " << totalC << endl;
        
        errors (nsteps * 0.5 * dh*dh);
    }
    
    void errors(double t)
    {
        double localLinf = 0., Linf = 0.;
        double localL1 = 0., L1 = 0.;
        double localL2 = 0., L2 = 0.;
        
        for (unsigned int iy=start_y; iy<end_y; iy++)
            for (unsigned int ix=start_x; ix<end_x; ix++)
            {
                // global position
                const int gix = rank_x * worksize_x + ix + k.getStart();
                const int giy = rank_y * worksize_y + iy + k.getStart();
                double error = cPing [ix + iy * workspace_x] - analytical ( (double)gix/(double)N, (double)giy/(double)N, t );
                localLinf = max ( localLinf, abs(error) );
                localL1 += abs (error);
                localL2 += error * error;
            }
        
        MPI_Reduce(&localLinf, &Linf, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&localL1,   &L1,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&localL2,   &L2,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if (rank==0)
        {
            L2 = sqrt(L2) / N;
            L1 /= N * N;
            cout << "Errors:\t" << Linf << "\t" << L1 << "\t" << L2 << endl;
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
    int Mx = 3;
    int My = 2;
    
    //=========================
    //      get rank information
    //      (processID and #processes)
    //=========================
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != Mx * My)
    {
      if (rank==0) cout << "Inconsistent number of MPI Processes\n";
      abort();
    }
    //printf("Hello from process %d of %d\n", rank, size);
    
    // Problem size
    int N = 120;
    
    // setup simulation
    DiffusionPSE <PSE2D> diffusion (rank, size, N, Mx, My);
    diffusion.setIC();
    
    // timer
    double timer = MPI_Wtime();
    
    // run simulation
    diffusion.simulate();
    
    // timer
    double runtime = MPI_Wtime() - timer;
    if (rank==0) cout << "Runtime: " << runtime << std::endl;
    
    //=========================
    //      cleanup MPI environment
    //=========================
    MPI_Finalize();
    
    return 0;
}
