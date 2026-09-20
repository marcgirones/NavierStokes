#include <iostream>
#include <cmath>


class Mesh 
{
private:
	int m_Nx {};
	int m_Ny {};
	double m_Lx {};
	double m_Ly {};
	double m_Lex {};
	double m_Ley {};

public:
	Mesh(int Nx, int Ny, double Lx, double Ly) : m_Nx { Nx }, m_Ny { Ny }, m_Lx { Lx }, m_Ly { Ly }, m_Lex { m_Lx/((double)m_Nx) }, m_Ley { m_Ly/((double)m_Ny) }
	{
		std::cout<<"mesh: elements ( "<<m_Nx<<", "<<m_Ny<<" ), dimensions of the element ( "<<m_Lex<<", "<<m_Ley<<" )\n"; 
	}
	void goeast (int i, int j, int* ie, int* je)
	{
		*je = j; 
		*ie = (i+1)%m_Nx;
	}
	double findxCenternode (int i)
	{
		return m_Lex/2 + m_Lex*i;
	}
	double findyCenternode (int j)
	{
		return m_Ley/2 + m_Ley*j;
	}
	double getLex ()
	{
		return m_Lex;
	}
	double getLey ()
	{
		return m_Ley;
	}
//		void nodeInfo (int N, int M)
//	{
//		std::cout<<" P("<<N<<","<<M<<") --> ("<<Mesh::findxPnode(N-1)<<", "<<Mesh::findyPnode(M-1)<<")\n U("<<N<<", "<<M<<") --> ("<<Mesh::findxUnode(N-1)<<", "<<Mesh::findyUnode(M-1)<<")\n V("<<N<<", "<<M<<") --> ("<<Mesh::findxVnode(N-1)<<", "<<Mesh::findyVnode(M-1)<<")\n";
//	}	
};
class Pmesh : public Mesh
{
public:
	Pmesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx,Ny,Lx,Ly)
	{
	}
	double findxPnode (int i)
	{
		return Mesh::findxCenternode(i);
	}
	double findyPnode (int j)
	{
		return Mesh::findyCenternode(j);
	}
};
class Umesh : public Mesh
{
public:
	Umesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx,Ny,Lx,Ly)
	{
	}
	double findxUnode (int i)
	{
		return (1+i)*Mesh::getLex();
	}
	double findyUnode (int j)
	{
		return Mesh::findyCenternode(j);
	}
};
class Vmesh : public Mesh
{
public:
	Vmesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx,Ny,Lx,Ly)
	{
	}
	double findxVnode (int i)
	{
		return Mesh::findxCenternode(i);
	}
	double findyVnode (int j)
	{
		return (1+j)*Mesh::getLey();
	}
};

int main ()
{
	int Nx = 10, Ny = 10; double Lx = 10.0, Ly = 10.0;
	Pmesh p { Nx, Ny, Lx, Ly};
	Umesh u { Nx, Ny, Lx, Ly};
	Vmesh v { Nx, Ny, Lx, Ly};
	std::cout<<"P(4,3) --> ("<<p.findxPnode(3)<<", "<<p.findyPnode(2)<<")\nU(4,3) --> ("<<u.findxUnode(3)<<", "<<u.findyUnode(2)<<")\nV(4,3) --> ("<<v.findxVnode(3)<<", "<<v.findyVnode(2)<<")\n";
	int j=0, ie=0, je=0;
	for (int i=0; i<Nx; i++)
	{
		p.goeast(i, j, &ie, &je);
		std::cout<<"Current node --> ("<<i<<", "<<j<<") ; east node --> ("<<ie<<", "<<je<<")\n";
	}
	return 0;
}

