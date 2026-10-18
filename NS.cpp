#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>     // atoi
#include <iomanip>     // setprecision
#include <string>      // to_string
#include "Mesh.h"
#include "NavierStokesTerms.h"
#include "Solver.h"
#define Lx 1.0
#define Ly 1.0



int main (int argc, char** argv)
{
	const int Nx = (argc>1) ? std::atoi(argv[1]) : 100; //default value N = 100
	const int Ny = Nx;

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
	
	/*td::cout<<"-------- Pressure initial field --------------------------------\n";
	p.print_field();
	std::cout<<"-------- Horizontal velocity initial field ----------------------\n";
	u.print_field();
	std::cout<<"-------- Vertical velocity initial field ------------------------\n";
	v.print_field();*/
	/**/	
	
	
	/*********************************************************
		PROVES INDEX, OPERATOR I TERMES DIFUSIU / CONVECTIU
	*********************************************************/

	/*int j=0, ie=0, je=0;
	double x=0, y =0, p_actual=0, p_east=0;
	double convU=0, convV=0, difU=0, difV=0, nu =1.51e-5;
	for (int i=0; i<Nx; i++)
	{
		//calcul de indexs tenint en compte halo, podem utilitzar el operator com veiem a baix
		ie = p.wrapX(i+1);
		je = p.wrapX(j);

		

		//pots asignar valors a un node [funcio operator]
		p(i+1,j) = 2*i;
		u(i+1,j) = 3*i;

		//ja no fa falta treballar amb index ie iw is.... simplement pots evaluar p(i+1,j)  [funcio operator]
		p_actual = p(i, j);
		p_east = p(i+1, j);

		//calcul termes convectius i difusius en u i v

		convU = FVM::computeConvectiveU(i, j, u, v);
		difU = FVM::computeDiffusiveU(i, j, u, nu);

		convV = FVM::computeConvectiveV(i, j, u, v);
		difV = FVM::computeDiffusiveV(i, j, u, v, nu);

		std::cout << "Convección U (" << i << ", " << j << "): " << convU << "\n";
		std::cout << "Difusión U   (" << i << ", " << j << "): " << difU << "\n";
		std::cout << "Convección V (" << i << ", " << j << "): " << convV << "\n";
		std::cout << "Difusión V   (" << i << ", " << j << "): " << difV << "\n";



		p.coord(i, j, &x, &y);
		std::cout<<"Current node --> ("<<i<<", "<<j<<") ; east node --> ("<<ie<<", "<<je<<")\n";
		//std::cout<<"Current node --> ("<<i<<", "<<j<<") ; coordinates --> ("<<x<<", "<<y<<")\n";
		std::cout<<"Current node value --> ("<<p_actual<<") ; East node value --> ("<<p_east<<")\n";
	}
		p.print_field();*/


	/*******************************************************
		CONVECTIVE AND DIFUSIVE TERM COMPUTATION
	********************************************************/
	double convU=0, convV=0, difU=0, difV=0, nu =1.51e-5;

	// Open CSV file for writing
	std::ofstream file("datos" + std::to_string(Nx) + ".csv");
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
	const double rho = 1.0;       // Set to the density required by the assignment
	const double dt = 1.0e-4;     // Example time step; check stability/convergence
	const int nSteps = 100;       // Number of time steps to compute

	Solver solver(Nx, Ny, Lx, Ly, nu, rho, dt);
	for (int n = 0; n < nSteps; ++n) {
		solver.step(u, v, p);
		if ((n + 1) % 10 == 0) {
			std::cout << "Completed step " << solver.getStep()
			          << ", t = " << solver.getTime() << "\n";
		}
	}

	// Final fields for inspection / verification.
	std::ofstream pressureFile("pressure_final.csv");
	std::ofstream uFile("u_final.csv");
	std::ofstream vFile("v_final.csv");
	pressureFile << std::setprecision(17);
	uFile << std::setprecision(17);
	vFile << std::setprecision(17);
	for (int j = 0; j < Ny; ++j) {
		for (int i = 0; i < Nx; ++i) {
			if (i) { pressureFile << ","; uFile << ","; vFile << ","; }
			pressureFile << p(i,j);
			uFile << u(i,j);
			vFile << v(i,j);
		}
		pressureFile << "\n"; uFile << "\n"; vFile << "\n";
	}

	// Spatial grid spacing
	const double dx = Lx / Nx;
	const double dy = Ly / Ny;

	// Verify incompressibility after the simulation
	double maxDiv = 0.0;
	double sumDiv2 = 0.0;

	for (int j = 0; j < Ny; ++j) {
		for (int i = 0; i < Nx; ++i) {

			const double div =
				(u(i, j) - u(i - 1, j)) / dx
			+ (v(i, j) - v(i, j - 1)) / dy;

			maxDiv = std::max(maxDiv, std::abs(div));
			sumDiv2 += div * div;
		}
	}

	const double rmsDiv =
		std::sqrt(sumDiv2 / (static_cast<double>(Nx) * Ny));

	std::cout << "\nIncompressibility verification:\n";
	std::cout << "Maximum divergence = " << maxDiv << '\n';
	std::cout << "RMS divergence     = " << rmsDiv << '\n';

	return 0;
}

