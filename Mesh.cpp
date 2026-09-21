#include "Mesh.h"
#include <iostream>
#include <vector>

Mesh::Mesh(int Nx, int Ny, double Lx, double Ly) : m_Nx { Nx }, m_Ny { Ny }, m_Lx { Lx }, m_Ly { Ly }, m_Lex { m_Lx/((double)m_Nx) }, m_Ley { m_Ly/((double)m_Ny) }, field(Nx*Ny)
{
	
	std::cout<<"mesh: elements ( "<<m_Nx<<", "<<m_Ny<<" ), dimensions of the element ( "<<m_Lex<<", "<<m_Ley<<" )\n"; 
	for( int ii=0; ii<field.size(); ii++){
		field[ii] = 0;
	}
	field[73] = 4;
}
void Mesh::goeast (int i, int j, int* ie, int* je)
{
	*je = j; 
	*ie = (i+1)%m_Nx;
}
double Mesh::findxCenternode (int i)
{
	return m_Lex/2 + m_Lex*i;
}
double Mesh::findyCenternode (int j)
{
	return m_Ley/2 + m_Ley*j;
}
double Mesh::findx (int i){return 0.0;}
double Mesh::findy (int j){return 0.0;}
void Mesh::coord (int i, int j, double* x, double* y){}
double Mesh::getLex ()
{
	return m_Lex;
}
double Mesh::getLey ()
{
	return m_Ley;
}
void Mesh::print_field ()
{
	int pos_vec=0; 
	for(int jj = m_Ny-1; jj>=0; jj--){
		for(int ii = 0; ii<m_Nx; ii++){
			pos_vec = jj*m_Nx + ii;
			std::cout<<field[pos_vec]<<" ";
		}
		//if(!pos_vec%m_Nx-1) std::cout<<"\n";
		std::cout<<"\n";
	}
}

Pmesh::Pmesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx,Ny,Lx,Ly){}
double Pmesh::findx (int i) 
{
	return Mesh::findxCenternode(i);
}
double Pmesh::findy (int j)
{
	return Mesh::findyCenternode(j);
}
void Pmesh::coord (int i, int j, double* x, double* y)
{
	*x = Pmesh::findx(i);
	*y = Pmesh::findy(j);
}

Umesh::Umesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx,Ny,Lx,Ly){}
double Umesh::findx (int i)
{
	return (1+i)*Mesh::getLex();
}
double Umesh::findy (int j)
{
	return Mesh::findyCenternode(j);
}
void Umesh::coord (int i, int j, double* x, double* y)
{
	*x = Umesh::findx(i);
	*y = Umesh::findy(j);
}


Vmesh::Vmesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx,Ny,Lx,Ly){}
double Vmesh::findx (int i)
{
	return Mesh::findxCenternode(i);
}
double Vmesh::findy (int j)
{
	return (1+j)*Mesh::getLey();
}
void Vmesh::coord (int i, int j, double* x, double* y)
{
	*x = Vmesh::findx(i);
	*y = Vmesh::findy(j);
}
