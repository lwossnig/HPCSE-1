#include <iostream>
#include <cassert>
#include <vector>
#include <fstream>
#include <cmath>
#include <complex>
//#define DEBUG
#define pi 3.14159265359

void update(std::vector<std::complex<double>> & , std::vector<std::complex<double>> & , const int , const double , const double );// update using the velocity and locality vector
void write_to_file(const char * , std::vector<std::complex<double>> & );

int main()
{
	int N = 1000;
	double h = 1/(double)N;
	std::cout << " Starting ... " << std::endl;
	std::vector<std::complex<double>> vortices(N);
	std::vector<std::complex<double>> velocities(N);
	
	std::cout << " Initialize particle positions... " << std::endl;
	for(int i=0;i<N;i++){
		vortices[i] = {-0.5+(i+0.5)*h,0};
		velocities[i] = {0,0};
#ifdef DEBUG
		std::cout << vortices[i] << std::endl;
#endif
	}
	double t = 0; const double dt = 0.0001; int no = 0, count = 0;
	char filename[160];
	/// updating the positions now
	while(t<3){
		if(count%100==0){
			sprintf(filename,"vortices%d.txt",no); // now same as "picture#nr"
			write_to_file(filename, vortices);
			no++;
		}
		update(vortices, velocities, N, h, dt);
#ifdef DEBUG
		for(int i=0;i<N;i++){
			std::cout << vortices[i] << std::endl;
		}
#endif
		count++;
		std::cout << "iteration with time t = " << t << std::endl;
		t += dt;

	}
	std::cout << " System in final state. t = " << t << std::endl;
	return 0;
}


void write_to_file(const char* filename, std::vector<std::complex<double>> &positions){
	std::ofstream out(filename);
	for(int i=0; i< positions.size()-1;i++){
		out << positions[i].real() << "    " << positions[i].imag() << std::endl;
	}
	out.close ();
}

void update(std::vector<std::complex<double>> &positions, std::vector<std::complex<double>> &velocities, const int N, const double h, const double dt)
{
	double Gamma_j(0); double x(0);
	std::complex<double> ic(0,1);
	for(int i=0; i<N; ++i){
		velocities[i] = {0,0};
		for(int j=0; j<N; ++j){
			if(j!=i){ 
				x = positions[j].real();
				Gamma_j = 4*x*h/std::sqrt(1-4*x*x);
				velocities[i] += (-1) * Gamma_j * ic / (2 * pi * ( positions[i] - positions[j] )); 
#ifdef DEBUG
				std::cout << velocities[i] << std::endl;
#endif
			}
		}
		positions[i] =  positions[i] + dt * std::conj(velocities[i]);
	}
}



