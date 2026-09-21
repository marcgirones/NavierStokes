#ifndef MESH_H
#define MESH_H

#include <vector>

class Mesh 
{
private:
	int m_Nx;
	int m_Ny;
	double m_Lx;
	double m_Ly;
	double m_Lex;
	double m_Ley;
	std::vector<double> field;
	

public:
	Mesh(int Nx, int Ny, double Lx, double Ly);
	void goeast (int i, int j, int* ie, int* je);
	double findxCenternode (int i);
	double findyCenternode (int j);
	virtual double findx (int i) = 0;
	virtual double findy (int j) = 0;
	virtual void coord (int i, int j, double* x, double* y) = 0;
	double getLex ();
	double getLey ();
	void print_field ();
};
class Pmesh : public Mesh
{
public:
	Pmesh(int Nx, int Ny, double Lx, double Ly);
	double findx (int i) override;
	double findy (int j) override;
	void coord (int i, int j, double* x, double* y) override;
};
class Umesh : public Mesh
{
public:
	Umesh(int Nx, int Ny, double Lx, double Ly);
	double findx (int i) override;
	double findy (int j) override;
	void coord (int i, int j, double* x, double* y) override;
};
class Vmesh : public Mesh
{
public:
	Vmesh(int Nx, int Ny, double Lx, double Ly);
	double findx (int i) override;
	double findy (int j) override;
	void coord (int i, int j, double* x, double* y) override;
};
#endif
