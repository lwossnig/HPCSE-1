/*
 *  main.cpp
 *
 *  Created by Christian Conti on 11/12/2012
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <sstream>

#include "Timer.h"
#include "ArgumentParser.h"
#include "Diffusion.h"


using namespace std;

int main(int argc, const char **argv)
{
	// setup helper objects
	ArgumentParser parser(argc, argv);
	Timer timer;
	
	// read input parameters from command line
	const int size = parser("-size").asInt(8);
    const Real dh = 1./(Real)(size-1);
	const Real D = parser("-D").asDouble(1.);
	Real tend = parser("-tend").asDouble(1.);

	// set additional helper constants
	const Real dt = parser("-dt").asDouble(dh*dh*.5/D);
	const int nsteps = parser("-nsteps").asInt(tend/dt);
    tend = nsteps*dt;
	
	// sets filenames for dumps of diagnostics
	stringstream euler_diag, adi_diag;
	euler_diag << "euler_diag_" << D << "_" << size << ".txt";
	adi_diag   << "adi_diag_"   << D << "_" << size << ".txt";
	
	
	cout << "Diffusion with ADI - " << nsteps << " steps\n";
	
	// create an object to compute diffusion with ADI and one for forward Euler
	// the latter is used to compare
	Diffusion diffusionEuler(size,D);
	Diffusion diffusionADI(size,D);
	diffusionEuler.dumpEuler("eulerIC.txt",0);
	diffusionADI.dumpADI("adiIC.txt",0);
	
	// compute diffusion with forward Euler
	timer.start();
	for (int t=0; t<nsteps; t++)
	{
		diffusionEuler.Euler(dt);
	}
	const double tEuler = timer.stop();
	
	// compute diffusion with ADI
	timer.start();
	for (int t=0; t<nsteps; t++)
	{
		diffusionADI.ADI(dt);
	}
	const double tADI = timer.stop();
	
	// output timings and grids
	cout << "Timings: " << D << "\t" << size << "\t" << tEuler << "\t" << tADI << endl;
	diffusionEuler.dumpEuler("euler.txt",tend);
	diffusionADI.dumpADI("adi.txt",tend);
	
	diffusionEuler.errors("eulerErr.txt", tend);
	diffusionADI.errors("adiErr.txt", tend);
	
	return 0;
} 
