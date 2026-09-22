#include <iostream>
#include <cmath>
#include "Mesh.h"
#define Nx 4
#define Ny 4
#define Lx 1.0
#define Ly 1.0



int main ()
{
	std::vector<double> init_field (Nx*Ny);
	for (int i=0; i<Nx*Ny; i++){init_field[i] = 0;}
	init_field[73] = 8; 
	Pmesh p { Nx, Ny, Lx, Ly};
	Umesh u { Nx, Ny, Lx, Ly};
	Vmesh v { Nx, Ny, Lx, Ly};	
	p.initialize([](double x, double y)
	{
		return 0;
	});

	
	u.initialize([](double x, double y)
	{
		double a = cos(2*M_PI*x)*sin(2*M_PI*y);
		if((a>0 && a>1e-10) || (a<0 && a<-1e-10)) return a;
		else return 0.0;
	});
	v.initialize([](double x, double y)
	{
		double a = -sin(2*M_PI*x)*cos(2*M_PI*y);
		if((a>0 && a>1e-10) || (a<0 && a<-1e-10)) return a;
		else return 0.0;
	});
	/*int j=0, ie=0, je=0;
	double x=0, y =0;
	for (int i=0; i<Nx; i++)
	{
		p.goeast(i, j, &ie, &je);
		p.coord(i, j, &x, &y);
		std::cout<<"Current node --> ("<<i<<", "<<j<<") ; east node --> ("<<ie<<", "<<je<<")\n";
		std::cout<<"Current node --> ("<<i<<", "<<j<<") ; coordinates --> ("<<x<<", "<<y<<")\n";
	}*/
	std::cout<<"-------- Pressure initial field --------------------------------\n";
	p.print_field();
	std::cout<<"-------- Horizontal velocity initial field ----------------------\n";
	u.print_field();
	std::cout<<"-------- Vertical velocity initial field ------------------------\n";
	v.print_field();
	return 0;
}

