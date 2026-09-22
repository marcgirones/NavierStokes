#include "Mesh.h"

Mesh::Mesh(int Nx, int Ny, double Lx, double Ly) 
    : m_Nx{Nx}, m_Ny{Ny}, m_Lx{Lx}, m_Ly{Ly}, 
      m_Lex{m_Lx / static_cast<double>(m_Nx)}, 
      m_Ley{m_Ly / static_cast<double>(m_Ny)}, 
      m_field(Nx * Ny, 0.0) 
{
    std::cout << "mesh: elements ( " << m_Nx << ", " << m_Ny 
              << " ), dimensions ( " << m_Lex << ", " << m_Ley << " )\n";
}

void Mesh::initialize (std::function<double(double,double)> init_cond)
{
	for(int ii = 0; ii<m_Nx; ii++){
		for(int jj = 0; jj<m_Ny; jj++){
			double x, y;
			this->coord(ii,jj,&x,&y);
			m_field[jj*m_Nx+ii] = init_cond(x,y);
		}
	}
}

double Mesh::findxCenternode(int i) const { return m_Lex / 2.0 + m_Lex * i; }
double Mesh::findyCenternode(int j) const { return m_Ley / 2.0 + m_Ley * j; }

void Mesh::print_field() const {
    for (int jj = m_Ny - 1; jj >= 0; --jj) {
        for (int ii = 0; ii < m_Nx; ++ii) {
            std::cout << (*this)(ii, jj) << " ";
        }
        std::cout << "\n";
    }
}

// Pmesh
Pmesh::Pmesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx, Ny, Lx, Ly) {}
double Pmesh::findx(int i) const { return findxCenternode(i); }
double Pmesh::findy(int j) const { return findyCenternode(j); }
void Pmesh::coord(int i, int j, double* x, double* y) const {
    *x = findx(i); *y = findy(j);
}

// Umesh
Umesh::Umesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx, Ny, Lx, Ly) {}
double Umesh::findx(int i) const { return (1 + i) * getLex(); }
double Umesh::findy(int j) const { return findyCenternode(j); }
void Umesh::coord(int i, int j, double* x, double* y) const {
    *x = findx(i); *y = findy(j);
}

// Vmesh
Vmesh::Vmesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx, Ny, Lx, Ly) {}
double Vmesh::findx(int i) const { return findxCenternode(i); }
double Vmesh::findy(int j) const { return (1 + j) * getLey(); }
void Vmesh::coord(int i, int j, double* x, double* y) const {
    *x = findx(i); *y = findy(j);
}