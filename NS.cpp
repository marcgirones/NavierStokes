#include <iostream>
#include <cmath>
#include "Mesh.h"
#define Nx 10
#define Ny 10
#define Lx 10.0
#define Ly 10.0



int main ()
{
	Pmesh p { Nx, Ny, Lx, Ly};
	Umesh u { Nx, Ny, Lx, Ly};
	Vmesh v { Nx, Ny, Lx, Ly};
	//std::cout<<"P(4,3) --> ("<<p.findx(3)<<", "<<p.findy(2)<<")\nU(4,3) --> ("<<u.findx(3)<<", "<<u.findy(2)<<")\nV(4,3) --> ("<<v.findx(3)<<", "<<v.findy(2)<<")\n";
	int j=0, ie=0, je=0;
	double x=0, y =0;
	for (int i=0; i<Nx; i++)
	{
		p.goeast(i, j, &ie, &je);
		p.coord(i, j, &x, &y);
		std::cout<<"Current node --> ("<<i<<", "<<j<<") ; east node --> ("<<ie<<", "<<je<<")\n";
		std::cout<<"Current node --> ("<<i<<", "<<j<<") ; coordinates --> ("<<x<<", "<<y<<")\n";
	}
		p.print_field();
	return 0;
}

