#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>     // atoi
#include <iomanip>     // setprecision
#include <string>      // to_string
#include "Mesh.h"
#include "NavierStokesTerms.h"
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

	// Close the file stream
	file.close();

	return 0;
}

