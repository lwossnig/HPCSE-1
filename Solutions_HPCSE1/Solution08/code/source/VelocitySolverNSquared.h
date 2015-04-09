
#ifndef __DirectVortexSolver__VelocitySolverNSquared__
#define __DirectVortexSolver__VelocitySolverNSquared__

#include "common.h"
#include "ArrayOfParticles.h"

class VelocitySolverNSquared
{
private:
    ArrayOfParticles & dstParticles;
    ArrayOfParticles & srcParticles; // this is read-only
    
    const int rank, size;
    
    const double wingSpan;
    
    
public:
	VelocitySolverNSquared(ArrayOfParticles & dstParticles, ArrayOfParticles & srcParticles, int rank, int size, const double wingSpan);
	~VelocitySolverNSquared(){}
    
    void ComputeVelocity();
    
    double timeC, timeT;
};


#endif /* defined(__DirectVortexSolver__VelocitySolverNSquared__) */
