/*
 *  Diffusion.h
 *
 *	Header file for Diffusion class
 *	Compute diffusion on a square domain [-1,1[^2
 *
 *	Parameters:
 *		Constructor - size: number of grid points per direction x and y
 *		Constructor - D: diffusion coefficient
 *		ADI - dt: timestep size
 *		Euler - dt: timestep size
 *
 *  Created by Christian Conti on 11/12/2012
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#include <string>

#include "Matrix2D.h"

typedef double Real;

class Diffusion
{
private:
	/*=================*/
	/* Private members */
	/*=================*/
	
	// Diffusion coefficient
	Real D;
	
	// Grid spacing and square of grid spacing
	Real dh, dh2;
	
	// Grid size
	unsigned int size;
	
	// Grids containing concentrations
	// grid_tmp, grid_tmp2 are necessary for ping-pong scheme
	Matrix2D<Real> grid, grid_tmp;
	
	
	/*=================*/
	/* Private methods */
	/*=================*/
	
	// Thomas algorithm for solving tridiagonal systems,
	//	- tailored to diffusion problem
	//	- dir==0: solve implicit diffusion step in x direction, dir==1: y direction
	//	- i: denotes row or column on which to compute the 1D implicit diffusion (depending on dir)
	//	- coeff: coefficient D*dt/2/(dh*dh)
	void Thomas(int dir, int i, Real coeff, Real dt);
    
    
	void _dump(std::string filename, double t);
	
public:
	/*=================*/
	/* Public methods  */
	/*=================*/
	
	// Constructor
	//	- requires grid size and diffusion coefficient D
	//	- sets initial conditions
	Diffusion(unsigned int size, Real D);
	
	// Advances the diffusion simulation by a single timestep dt using ADI
	//	assumes dirichlet boundary conditions of 0
	void ADI(Real dt);
	
	// Advances the diffusion simulation by a single timestep dt using forward Euler
	//	assumes dirichlet boundary conditions of 0
	void Euler(Real dt);
	
	
	/*=================*/
	/* Utilities       */
	/*=================*/
	
	// saves the content of grid to file filename
	void dumpEuler(std::string filename, double t);
	void dumpADI(std::string filename, double t);
	
	// saves the errors with respect to analytical solution
	void errors(std::string filename, double t);
    
    double analytical(unsigned int ix, unsigned int iy, double t);
};
