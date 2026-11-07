#include "Mesh.hpp"

Mesh::Mesh(int Nx, int Ny, double Lx, double Ly) 
    : m_Nx{Nx}, m_Ny{Ny}, m_Lx{Lx}, m_Ly{Ly}, 
      m_Lex{m_Lx / static_cast<double>(m_Nx)}, 
      m_Ley{m_Ly / static_cast<double>(m_Ny)}, 
      m_field(Nx * Ny, 0.0) 
{
}
///////////////////////////////////////////////////////////////////////////////////////////////
///Initialitzation of the field of the object /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

void Mesh::initialize (std::function<double(double,double)> init_cond)
{
	//initialize fills the values inside any field inside the Mesh class (Nx*Ny vector) using a function regarding the position of the node
	// where field[i][j] = f(x,y) being x and y the physical position of the center of the node (i,j)

	for(int ii = 0; ii<m_Nx; ii++){
		for(int jj = 0; jj<m_Ny; jj++){
			double x, y;
			this->coord(ii,jj,&x,&y);
			m_field[jj*m_Nx+ii] = init_cond(x,y);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////
///Class functions ////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

double Mesh::findxCenternode(int i) const { return m_Lex / 2.0 + m_Lex * i; } //finds the non-staggered center of the node i in the x axis
double Mesh::findyCenternode(int j) const { return m_Ley / 2.0 + m_Ley * j; } //finds the non-staggered center of the node i in the y axis

// Set of functions that recive the position on the matrix (i,j) or the vector (i) and returns the index of the desired neighbour node 
// east, west, north or south in the vector configuration

int Mesh::goeast (int i, int j) const{
	return j*m_Nx + (i+1+m_Nx)%m_Nx;
}
int Mesh::goeast (int i) const{
	return (i/m_Nx)*m_Nx + (i+1)%m_Nx;
}

int Mesh::gowest (int i, int j) const{
	return j*m_Nx + (i-1+m_Nx)%m_Nx;
}
int Mesh::gowest (int i) const{
	return (i/m_Nx)*m_Nx + (i-1+m_Nx)%m_Nx;
}
int Mesh::gonorth (int i, int j) const{
	return ((j+1)%m_Ny)*m_Nx + i;
}
int Mesh::gonorth (int i) const{
	return i%m_Nx + (((i/m_Nx)+1)%m_Ny)*m_Nx;
}
int Mesh::gosouth (int i, int j) const{
	return ((j-1+m_Ny)%m_Ny)*m_Nx + i;
}
int Mesh::gosouth (int i) const{
	return i%m_Nx + (((i/m_Nx)-1+m_Ny)%m_Ny)*m_Nx; 
}


double Mesh::max_val() {
    //Return the maximum value in the field of this object

    double maxval = 0.0;

    for (int i = 0; i < m_Nx * m_Ny; i++) {
        const double absval = std::abs(m_field[i]);
        absval > maxval ? maxval = absval : maxval = maxval;
    }

    return maxval;
}

void Mesh::print_field() const {
    //Prints the values of the field in a physical order, meaning that the value printed in the bottom left is the (0) or (0,0) of the field

    for (int jj = m_Ny - 1; jj >= 0; --jj) {
        for (int ii = 0; ii < m_Nx; ++ii) {
            std::cout << (*this)(ii, jj) << " ";
        }
        std::cout << "\n";
    }
}

double Mesh::getFieldValue (int i) const{
	return m_field[i];
}
void Mesh::setFieldValue (int i, double val){
	m_field[i] = val;
}

///////////////////////////////////////////////////////////////////////////////////////////////
///Pmesh //////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

//Subclass of the Mesh class that is centered to the phisical discretitzation of the domain, as it would be in a non-staggered mesh

Pmesh::Pmesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx, Ny, Lx, Ly) {}
double Pmesh::findx(int i) const { return findxCenternode(i); }
double Pmesh::findy(int j) const { return findyCenternode(j); }
void Pmesh::coord(int i, int j, double* x, double* y) const {
    *x = findx(i); *y = findy(j);
}

///////////////////////////////////////////////////////////////////////////////////////////////
///Umesh //////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

//Subclass of the Mesh class, which field is moved 0.5 dx in the x axis for the staggered horitzontal velocity 

Umesh::Umesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx, Ny, Lx, Ly) {}
double Umesh::findx(int i) const { return (1 + i) * getLex(); }
double Umesh::findy(int j) const { return findyCenternode(j); }
void Umesh::coord(int i, int j, double* x, double* y) const {
    //Return the physical position of the node in the staggered node

    *x = findx(i); *y = findy(j);
}
double Umesh::gradi (int ii) const{
	//Returns the local gradient at the node (i)
	
	return (getFieldValue(ii) - getFieldValue(gowest(ii)))/getLex();
}

///////////////////////////////////////////////////////////////////////////////////////////////
///Umesh //////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

//Subclass of the Mesh class, which field is moved 0.5 dy in the y axis for the staggered vertical velocity

Vmesh::Vmesh(int Nx, int Ny, double Lx, double Ly) : Mesh(Nx, Ny, Lx, Ly) {}
double Vmesh::findx(int i) const { return findxCenternode(i); }
double Vmesh::findy(int j) const { return (1 + j) * getLey(); }
void Vmesh::coord(int i, int j, double* x, double* y) const {
    //Return the physical position of the node in the staggered node

    *x = findx(i); *y = findy(j);
}
double Vmesh::gradi (int ii) const{
	//Returns the local gradient at the node (i)
	
	return (getFieldValue(ii) - getFieldValue(gosouth(ii)))/getLey();
}
