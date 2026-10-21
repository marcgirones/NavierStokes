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
		PART A CONVECTIVE AND DIFUSIVE TERM COMPUTATION
	********************************************************/
	double convU=0, convV=0, difU=0, difV=0;

	// Open CSV file for writing
	std::ofstream file("datos" + std::to_string(Nx) + ".csv");
	if (!file) {
		std::cerr << "Error: could not open Part A CSV output.\n";
		return 1;
	}
	file << std::setprecision(17);

	// Write header row
	file << "convU,convV,difU,difV\n";

	for (int i=0; i<Nx; i++)
	{
		for (int j=0; j<Ny; j++)
		{
			// Compute convective and diffusive terms
			convU = FVM::computeConvectiveU(i, j, u, v);
			difU = FVM::computeDiffusiveU(i, j, u, nu);

			convV = FVM::computeConvectiveV(i, j, u, v);
			difV = FVM::computeDiffusiveV(i, j, u, v, nu);

			// Write row values separated by commas
			file << convU << "," << convV << "," << difU << "," << difV << "\n";

		}
	}

	// Close the Part A verification file. The original verification loop above
	// is intentionally preserved unchanged.
	file.close();
}
