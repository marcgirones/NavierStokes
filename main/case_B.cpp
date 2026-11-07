#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>     // atoi
#include <iomanip>     // setprecision
#include <string>      // to_string
#include <algorithm>   // max
#include "Mesh.hpp"
#include "NavierStokesTerms.hpp"
#include "Solver.hpp"
#define Lx 1.0
#define Ly 1.0



int main (int argc, char** argv)
{
	const int Nx = (argc>1) ? std::atoi(argv[1]) : 100; //default value N = 100
	const int Ny = Nx;

	// Simulation parameters: N, final time, viscosity, time step
	const double tFinal = (argc > 2) ? std::atof(argv[2]) : 0.01;
	const double nu = (argc > 3) ? std::atof(argv[3]) : 1.51e-5;
	const double rho = 1.0;

	if (Nx < 1 || tFinal <= 0.0 || nu < 0.0) {
		std::cerr << "Error: N and dt/tFinal must be positive, and nu must be non-negative.\n";
		return 1;
	}

	Pmesh p { Nx, Ny, Lx, Ly};
	Umesh u { Nx, Ny, Lx, Ly};
	Vmesh v { Nx, Ny, Lx, Ly};

	/************************************
		Initialize the field
	************************************/

	p.initialize([](double x, double y)
	{
		return 0;
	});


	u.initialize([](double x, double y)
	{
		//double a = cos(2*M_PI*x)*sin(2*M_PI*y);
		//double a = cos(x);
		return 0.0;
	});
	v.initialize([](double x, double y)
	{
		//double a = -sin(2*M_PI*x)*cos(2*M_PI*y);
		//double a = sin(y);
		return 0.0;
	});

	u.setFieldValue(Nx*Ny/2+50, 100*1.0/Nx);

	/*******************************************************
		PART B: PRESSURE-VELOCITY COUPLING AND TIME LOOP
	*******************************************************/
	std::ofstream file1("results/case_B_results/upFieldDivergence.csv");
	std::ofstream file2("results/case_B_results/unFieldDivergence.csv");


	Solver solver(Nx, Ny, Lx, Ly, nu, rho);

	double maxdiv = 0, totaldiv=0;

	for(int i=0; i<Nx*Ny; i++){

		//Check the initial local and global divergence of the velocity field 

		double divu = u.gradi(i);
		double divv = v.gradi(i);
		double div = std::abs(divu+divv);
		totaldiv+=divu+divv;
		if(div>maxdiv) maxdiv = div;
		if(i%Nx==Nx-1) file1<<divu+divv<<"\n";
		else file1<<divu+divv<<",";
	}
	std::cout<<"maxdiv (up) "<<maxdiv<<" Total divergence of the field "<<totaldiv<<"\n";

	//Solve the poisson probelm

	solver.solvePoissonFFT(u, v, p);


	maxdiv = 0; totaldiv = 0;
	for(int i=0; i<Nx*Ny; i++){
		
		//Check the local and global divergence of the u^(n+1) velocity field

		double divu = u.gradi(i);
		double divv = v.gradi(i);
		double div = std::abs(divu+divv);
		totaldiv+= divu+divv;
		if(div>maxdiv) maxdiv = div;
		if(i%Nx==Nx-1) file2<<divu+divv<<"\n";
		else file2<<divu+divv<<",";

	}
	std::cout<<"maxdiv (u^n+1) "<<maxdiv<<" Total divergence of the field "<<totaldiv<<"\n";
	file1.close(); file2.close();
}
