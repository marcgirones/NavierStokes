#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>     // atoi
#include <iomanip>     // setprecision
#include <string>      // to_string
#include <algorithm>   // max
#include <vector>
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
	// Re = U*L/nu = 100 with U = 1 (max |u| of the Taylor-Green vortex) and L = 1
	const double nu = (argc > 3) ? std::atof(argv[3]) : 0.01;
	// Save a snapshot every 'saveEvery' steps (1 = every time step).
	// The initial state and the final state are always saved.
	const int saveEvery = (argc > 4) ? std::max(1, std::atoi(argv[4])) : 1;
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
	});
	v.initialize([](double x, double y)
	{
		double a = -sin(2*M_PI*x)*cos(2*M_PI*y);
		return a;
	});



	/*******************************************************
		PART A CONVECTIVE AND DIFUSIVE TERM COMPUTATION
	********************************************************/
	double convU=0, convV=0, difU=0, difV=0;

	for (int i=0; i<Nx*Ny; i++)
	{
		// Compute convective and diffusive terms
		convU = FVM::computeConvectiveU(i, u, v);
		difU = FVM::computeDiffusiveU(i, u, nu);

		convV = FVM::computeConvectiveV(i, u, v);
		difV = FVM::computeDiffusiveV(i, u, v, nu);

	}

	/*******************************************************
		PART B/C: PRESSURE-VELOCITY COUPLING AND TIME LOOP
	*******************************************************/
	Solver solver(Nx, Ny, Lx, Ly, nu, rho);

	/*******************************************************
		HISTORY OUTPUT (all states, one block per saved step)
		history_u.bin / history_v.bin / history_p.bin:
		  raw float64, each snapshot = Ny*Nx values, row-major
		  (index j*Nx+i), snapshots stored consecutively.
		history_times.csv: step,t,dt for every snapshot.
	********************************************************/
	std::ofstream hU("history_u.bin", std::ios::binary);
	std::ofstream hV("history_v.bin", std::ios::binary);
	std::ofstream hP("history_p.bin", std::ios::binary);
	std::ofstream hT("history_times.csv");
	if (!hU || !hV || !hP || !hT) {
		std::cerr << "Error: could not open history output files.\n";
		return 1;
	}
	hT << std::setprecision(17) << "step,t,dt\n";
	
	int end =0;
	
	//generation of the binary files to store all transitions

	while (!end) {
 
		if (solver.getStep() % saveEvery == 0 || solver.getTime() >= tFinal) {
			for (int k = 0; k < Nx*Ny; k++) {
				const double uk = u.getFieldValue(k);
				const double vk = v.getFieldValue(k);
				const double pk = p.getFieldValue(k);
				hU.write(reinterpret_cast<const char*>(&uk), sizeof(double));
				hV.write(reinterpret_cast<const char*>(&vk), sizeof(double));
				hP.write(reinterpret_cast<const char*>(&pk), sizeof(double));
			}
			hT << solver.getStep() << "," << solver.getTime() << ","<< solver.getdt() << "\n";
		}
 
		if (solver.getTime() >= tFinal) end=1;
 		else{
			const double dtMax = tFinal - solver.getTime();
			solver.step(u, v, p, dtMax);
 
			std::cout << "Completed step " << solver.getStep()<< ", dt = " << solver.getdt()<< ", t = " << solver.getTime() << "\n";
		}
	}

	hU.close(); hV.close(); hP.close(); hT.close();

	//Generation of a .csv for the last snapshot
	std::ofstream pressureNFile("pressure_final_" + std::to_string(Nx) + ".csv");
	std::ofstream uNFile("u_final_" + std::to_string(Nx) + ".csv");
	std::ofstream vNFile("v_final_" + std::to_string(Nx) + ".csv");

	if (!pressureNFile || !uNFile || !vNFile) {
		std::cerr << "Error: could not open one or more final-field CSV files.\n";
		return 1;
	}

	pressureNFile << std::setprecision(17);
	uNFile << std::setprecision(17);
	vNFile << std::setprecision(17);

	for (int j = 0; j < Ny; ++j) {
		for (int i = 0; i < Nx; ++i) {
			if (i) {
				pressureNFile << ","; uNFile << ","; vNFile << ",";
			}
			pressureNFile << p.getFieldValue(i,j);
			uNFile << u.getFieldValue(i,j);
			vNFile << v.getFieldValue(i,j);
		}
		pressureNFile << "\n"; uNFile << "\n"; vNFile << "\n";
	}

	pressureNFile.close();
	uNFile.close();
	vNFile.close();

	return 0;
}
