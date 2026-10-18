#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>     // atoi
#include <iomanip>     // setprecision
#include <string>      // to_string
#include <algorithm>   // max
#include "Mesh.h"
#include "NavierStokesTerms.h"
#include "Solver.h"
#define Lx 1.0
#define Ly 1.0



int main (int argc, char** argv)
{
	const int Nx = (argc>1) ? std::atoi(argv[1]) : 100; //default value N = 100
	const int Ny = Nx;

	// Simulation parameters: N, final time, viscosity, time step
	const double tFinal = (argc > 2) ? std::atof(argv[2]) : 0.01;
	const double nu = (argc > 3) ? std::atof(argv[3]) : 1.51e-5;
	const double dt = (argc > 4) ? std::atof(argv[4]) : 1.0e-4;
	const double rho = 1.0;

	if (Nx < 1 || tFinal <= 0.0 || nu < 0.0 || dt <= 0.0) {
		std::cerr << "Error: N and dt/tFinal must be positive, and nu must be non-negative.\n";
		return 1;
	}

	const int nSteps = static_cast<int>(std::round(tFinal / dt));
	if (nSteps < 1 || std::abs(nSteps * dt - tFinal) > 1e-12) {
		std::cerr << "Error: tFinal must be a positive integer multiple of dt.\n";
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

	/*******************************************************
		PART B/C: PRESSURE-VELOCITY COUPLING AND TIME LOOP
	*******************************************************/
	Solver solver(Nx, Ny, Lx, Ly, nu, rho, dt);
	for (int n = 0; n < nSteps; ++n) {
		solver.step(u, v, p);
		if ((n + 1) % 10 == 0) {
			std::cout << "Completed step " << solver.getStep()
			          << ", t = " << solver.getTime() << "\n";
		}
	}

	// Final fields for inspection / verification.
	// Original filenames are preserved.
	std::ofstream pressureFile("pressure_final.csv");
	std::ofstream uFile("u_final.csv");
	std::ofstream vFile("v_final.csv");

	// Additional N-specific files prevent results from different meshes
	// overwriting each other during a convergence study.
	std::ofstream pressureNFile("pressure_final_" + std::to_string(Nx) + ".csv");
	std::ofstream uNFile("u_final_" + std::to_string(Nx) + ".csv");
	std::ofstream vNFile("v_final_" + std::to_string(Nx) + ".csv");

	if (!pressureFile || !uFile || !vFile ||
	    !pressureNFile || !uNFile || !vNFile) {
		std::cerr << "Error: could not open one or more final-field CSV files.\n";
		return 1;
	}

	pressureFile << std::setprecision(17);
	uFile << std::setprecision(17);
	vFile << std::setprecision(17);
	pressureNFile << std::setprecision(17);
	uNFile << std::setprecision(17);
	vNFile << std::setprecision(17);

	for (int j = 0; j < Ny; ++j) {
		for (int i = 0; i < Nx; ++i) {
			if (i) {
				pressureFile << ","; uFile << ","; vFile << ",";
				pressureNFile << ","; uNFile << ","; vNFile << ",";
			}
			pressureFile << p(i,j);
			uFile << u(i,j);
			vFile << v(i,j);

			pressureNFile << p(i,j);
			uNFile << u(i,j);
			vNFile << v(i,j);
		}
		pressureFile << "\n"; uFile << "\n"; vFile << "\n";
		pressureNFile << "\n"; uNFile << "\n"; vNFile << "\n";
	}

	pressureFile.close();
	uFile.close();
	vFile.close();
	pressureNFile.close();
	uNFile.close();
	vNFile.close();

	return 0;
}
