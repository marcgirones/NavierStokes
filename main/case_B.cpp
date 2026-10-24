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
		PROVES SET FIELD
	************************************/

	p.initialize([](double x, double y)
	{
		return 0;
	});


	u.initialize([](double x, double y)
	{
		double a = cos(2*M_PI*x)*sin(2*M_PI*y);
		return a;
		//if((a>0 && a>1e-10) || (a<0 && a<-1e-10)) return a;
		//else return 0.0;
	});
	v.initialize([](double x, double y)
	{
		double a = -sin(2*M_PI*x)*cos(2*M_PI*y);
		return a;
		//if((a>0 && a>1e-10) || (a<0 && a<-1e-10)) return a;
		//else return 0.0;
	});



	/*******************************************************
		PART B/C: PRESSURE-VELOCITY COUPLING AND TIME LOOP
	*******************************************************/
	Solver solver(Nx, Ny, Lx, Ly, nu, rho);
		
	solver.step(u, v, p, 0.1);

	double maxdiv=0;

	for(int ii = 0; ii<Nx*Ny; ii++){
		dudx = (u.field(ii) + u.field(goeast(ii))/(Nx/Lx);
		dvdy = (v.field(ii) + v.field(gonorth(ii))/(Ny/Ly);
	}
}
