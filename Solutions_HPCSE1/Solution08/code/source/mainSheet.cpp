#include <fstream>
#include "common.h"
#include "ArrayOfParticles.h"
#include "VelocitySolverNSquared.h"

#define WEAK_SCALING

void WriteToFile(ArrayOfParticles & allParticles, const int Np, const char * sFileName)
{
    FILE * f = fopen(sFileName, "w");
    
    fprintf(f,"%s,%s,%s,%s,%s,%s,%s\n","x","y","z","vx","vy","vz","w");
    for(int i=0; i<Np;i++)
        fprintf(f,"%f,%f,%f,%f,%f,%f,%f\n", allParticles.x[i],allParticles.y[i],0.0,
                allParticles.u[i],allParticles.v[i],0.0,
                allParticles.gamma[i]);
    fclose(f);
}

void WriteTimeToFile(const int Nranks, const double time, const char * sFileName)
{
    FILE * f = fopen(sFileName, "a");
    
    fprintf(f,"%d\t%f\n", Nranks,time);
    
    fclose(f);
}
						 
int main (int argc, char ** argv)
{
    MPI_Init(&argc,&argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // number of processes
    if (rank==0)
        std::cout << "Running with " << size << " MPI processes\n";
    
    // timer setup
    Timer timerIC, timerSim;
	
	// time integration setup
	const double dt = 0.001;
    const double tfinal = 0.01;//2.5;
	const int itmax = 1e9;
	const int ndump = std::floor(0.01/dt);
    const bool bAnimation = false;
    const bool bVerbose = false;
	
    timerIC.start();
	// initialization
#ifndef WEAK_SCALING
    //size_t Np = (10000/size)*size; // ensures that the number of particles is divisible by the number of workers
    size_t Np = (100000/size)*size; // ensures that the number of particles is divisible by the number of workers
#else
    //size_t Np = 10000*sqrt(size);
    size_t Np = 100000*sqrt(size);
#endif
    
    // each worker gets part of the whole array
    size_t NpProcess = Np/size;
    ArrayOfParticles dstParticles(NpProcess), srcParticles(NpProcess);
    ArrayOfParticles allParticles(rank==0 ? Np : 0);
    
    const double wingSpan = 10.;
	const double dx = wingSpan/Np;
    double totGamma=0;
    
	for(size_t i=0; i<NpProcess; i++)
	{
        // assumptions
        // b = wingSpan
        // Gamma_s = max_circ
        const double max_circ = 1.;
        
		double x = -wingSpan*.5 + ((double)i+(double)rank*(double)NpProcess+.5)*dx;
		dstParticles.x[i] = x;
		dstParticles.y[i] = 0.0;
        double g = max_circ*(x/(wingSpan*.5))/sqrt(1.-x*x/((wingSpan*.5)*(wingSpan*.5)));
		dstParticles.gamma[i] = g*dx;
		totGamma += dstParticles.gamma[i];
        
        // initially, each rank has to compute the interactions with its own particles
		srcParticles.x[i] = x;
		srcParticles.y[i] = 0.0;
		srcParticles.gamma[i] = g*dx;
	}
    MPI_Reduce(rank==0 ? MPI_IN_PLACE : &totGamma, &totGamma, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    timerIC.stop();
    
    if (rank==0)
    {
        std::cout << "Number of particles: " << Np << std::endl;
        std::cout << "Number of particles per process: " << NpProcess << std::endl;
        std::cout << "Initial circulation: " << totGamma << std::endl;
        std::cout << "IC time: " << timerIC.get_timing() << std::endl;
    }
    
    
    VelocitySolverNSquared VelocitySolver(dstParticles, srcParticles, rank, size, wingSpan);
	
    timerSim.start();
    double t=0;
	for (int i=0; i<itmax; i++)
	{
		if(i%ndump==0 && rank==0 && bVerbose)
			std::cout << "Iteration " << i << " time " << t << std::endl;

		// compute velocities
		dstParticles.ClearVelocities();
        
		VelocitySolver.ComputeVelocity();
        
        // dump the particles
        if(i%ndump==0 && bAnimation)
        {
            // need to gather all particles here
            MPI_Gather(dstParticles.x    , NpProcess, MPI_DOUBLE, &allParticles.x[rank*NpProcess]    , NpProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(dstParticles.y    , NpProcess, MPI_DOUBLE, &allParticles.y[rank*NpProcess]    , NpProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(dstParticles.u    , NpProcess, MPI_DOUBLE, &allParticles.u[rank*NpProcess]    , NpProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(dstParticles.v    , NpProcess, MPI_DOUBLE, &allParticles.v[rank*NpProcess]    , NpProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(dstParticles.gamma, NpProcess, MPI_DOUBLE, &allParticles.gamma[rank*NpProcess], NpProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            
            if (rank==0)
            {
                char buf[500];
                sprintf(buf, "vortexsheet_%5.5d.csv", i);
                WriteToFile(allParticles,Np,buf);
            }
        }
		
		// advect Euler
		dstParticles.AdvectEuler(dt);
        srcParticles = dstParticles;
		
        if (rank==0)
            t += dt;
        
        MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		if (t>tfinal)
        {
            std::cout << "Bye from rank " << rank << std::endl;
            break;
        }
    }
    timerSim.stop();
    
    if (rank==0)
    {
        char buf[500];
        sprintf(buf, "timing.dat");
        WriteTimeToFile(size,timerSim.get_timing(),buf);
        std::cout << "#Ranks, Time - " << size << "\t" << timerSim.get_timing() << "\t( " << VelocitySolver.timeT << "\t" << VelocitySolver.timeC << " )\n";
    }
    
    MPI_Finalize();

	return 0;
}


