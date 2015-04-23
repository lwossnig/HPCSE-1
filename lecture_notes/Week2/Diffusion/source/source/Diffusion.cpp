/*
 *  Diffusion.cpp
 *
 *  Created by Christian Conti on 11/12/2012
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <vector>

#include "Diffusion.h"

using namespace std;

// Thomas algorithm for solving tridiagonal systems
void Diffusion::Thomas(int dir, int i, Real coeff, Real dt)
{
	// solves tridiagonal system of type Ax=v
	
#ifndef NDEBUG
	// ensure that we are only treating x and y directions
	assert(dir==0 || dir==1);
#endif
	
	// temporary vector to store the modified main diagonal of the tridiagonal matrix
    
    // PRECOMPUTE ALL COEFFICIENTS - THEY ARE THE SAME ACROSS TIME AND SPACE!!!
	vector<Real> b_tmp(size, 2.*(1./dt + coeff));
    Real a = -coeff;
    Real c = -coeff;
	
	if (dir==0)
	{
		// x direction
		//	grid_tmp contains solution at time t+dt/2
		//	grid contains solution at time t+dt
		
		// compute modified coefficients
		for (unsigned int ix=2; ix<size-1; ix++)
		{
			const Real m = a / b_tmp[ix-1];
			b_tmp[ix] -= m * c;
			grid_tmp(ix,i) -= m * grid_tmp.Read(ix-1,i);
		}
		
		// back substitution phase
		grid(size-2,i) = grid_tmp.Read(size-2,i) / b_tmp[size-2];
		
		for (int ix=size-3; ix>=1; ix--)
			grid(ix,i) = (grid_tmp.Read(ix,i) - c * grid.Read(ix+1,i)) / b_tmp[ix];
	}
	else
	{
		// y direction
		//	grid_tmp contains solution at time t+dt/2
		//	grid contains solution at time t+dt
		
		// compute modified coefficients
		for (unsigned int iy=2; iy<size-1; iy++)
		{
			const Real m = a / b_tmp[iy-1];
			b_tmp[iy] -= m * c;
			grid_tmp(i,iy) -= m * grid_tmp.Read(i,iy-1);
		}
		
		// back substitution phase
		grid(i,size-2) = grid_tmp.Read(i,size-2) / b_tmp[size-2];
		
		for (int iy=size-3; iy>=1; iy--)
			grid(i,iy) = (grid_tmp.Read(i,iy) - c * grid.Read(i,iy+1)) / b_tmp[iy];
	}
}


Diffusion::Diffusion(unsigned int size, Real D) : grid(size,size), grid_tmp(size,size), size(size), D(D), dh(1./(Real)(size-1))
{
	dh2 = dh*dh;
	
	// set initial conditions
	// parallel for, for NUMA 1st touch policy
#pragma omp parallel for
	for (unsigned int iy=0; iy<size; iy++)
		for (unsigned int ix=0; ix<size; ix++)
		{
            grid(ix,iy)      = analytical(ix,iy,0);
            grid_tmp(ix,iy)  = analytical(ix,iy,0);
		}
}

// Compute a single diffusion step with ADI
void Diffusion::ADI(Real dt)
{
	const Real coeff = D/dh2;
    
	// Explicit y
#pragma omp parallel for
	for (unsigned int iy=1; iy<size-1; iy++)
		for (unsigned int ix=1; ix<size-1; ix++)
			grid_tmp(ix,iy) = coeff * grid.Read(ix,iy-1) + 2. * (1./dt - coeff) * grid.Read(ix,iy) + coeff * grid.Read(ix,iy+1);
	
	// Implicit x
#pragma omp parallel for
	for (unsigned int iy=1; iy<size-1; iy++)
		Thomas(0,iy,coeff,dt);
    
	// Explicit x
#pragma omp parallel for
	for (unsigned int iy=1; iy<size-1; iy++)
		for (unsigned int ix=1; ix<size-1; ix++)
			grid_tmp(ix,iy) = coeff * grid.Read(ix-1,iy) + 2. * (1./dt - coeff) * grid.Read(ix,iy) + coeff * grid.Read(ix+1,iy);
	
	// Implicit y
#pragma omp parallel for
	for (unsigned int ix=1; ix<size-1; ix++)
		Thomas(1,ix,coeff,dt);
}

// Compute a single diffusion step with forward (explicit) Euler
void Diffusion::Euler(Real dt)
{
	const Real coeff = D*dt/dh2;
	
#pragma omp parallel for
	for (unsigned int iy=1; iy<size-1; iy++)
		for (unsigned int ix=1; ix<size-1; ix++)
			grid_tmp(ix,iy) = grid.Read(ix,iy) + coeff * (grid.Read(ix-1,iy) + grid.Read(ix,iy-1) - 4 * grid.Read(ix,iy) + grid.Read(ix+1,iy) + grid.Read(ix,iy+1));
	
	grid.swap(grid_tmp);
}

void Diffusion::_dump(string filename, double t)
{
	return;
	ofstream myfile(filename);
	if (myfile.is_open())
	{
		for (unsigned int iy=0; iy<size; iy++)
		{
			for (unsigned int ix=0; ix<size; ix++)
				myfile << fixed << setprecision(8) << grid(ix,iy) << "\t";
			myfile << endl;
		}
	}
   	myfile.close();
}

void Diffusion::errors(string filename, double t)
{
    Real Linf = 0.;
    Real L1 = 0.;
    Real L2 = 0.;
	
	ofstream myfile(filename, fstream::app);
	
	for (unsigned int iy=0; iy<size; iy++)
		for (unsigned int ix=0; ix<size; ix++)
        {
            if (ix==0 || iy==0 || ix==size-1 || iy==size-1)
                assert(grid(ix,iy)==0);
            Real error = grid(ix,iy)-analytical(ix,iy,t);
            Linf = max(Linf,abs(error));
            L1 += abs(error);
            L2 += error*error;
        }
    
    L2 = sqrt(L2)/size;
	L1 /= size*size;
    myfile << Linf << " " << L1 << " " << L2 << endl;
}

void Diffusion::dumpEuler(string filename, double t)
{
    _dump(filename,t);
}

void Diffusion::dumpADI(string filename, double t)
{
    _dump(filename,t);
}

double Diffusion::analytical(unsigned int ix, unsigned int iy, double t)
{
    return sin((double)ix*dh*2.*M_PI)*sin((double)iy*dh*2.*M_PI)*exp(-8.*D*M_PI*M_PI*t);
}
