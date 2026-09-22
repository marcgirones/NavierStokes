#ifndef MESH_H
#define MESH_H

#include <vector>
#include <iostream>

class Mesh {
protected:
    int m_Nx;
    int m_Ny;
    double m_Lx;
    double m_Ly;
    double m_Lex;
    double m_Ley;
    std::vector<double> m_field;

    
public:
    Mesh(int Nx, int Ny, double Lx, double Ly);
    virtual ~Mesh() = default;

	// Ajuste de índice periódico para 1D
    int wrapX(int i) const { return (i % m_Nx + m_Nx) % m_Nx; }
    int wrapY(int j) const { return (j % m_Ny + m_Ny) % m_Ny; }


    // Acceso periódico a los datos del campo
    double operator()(int i, int j) const {
        return m_field[wrapY(j) * m_Nx + wrapX(i)];
    }

    double& operator()(int i, int j) {
        return m_field[wrapY(j) * m_Nx + wrapX(i)];
    }

    // Getters de dimensiones
    int getNx() const { return m_Nx; }
    int getNy() const { return m_Ny; }
    double getLex() const { return m_Lex; }
    double getLey() const { return m_Ley; }

    double findxCenternode(int i) const;
    double findyCenternode(int j) const;

    virtual double findx(int i) const = 0;
    virtual double findy(int j) const = 0;
    virtual void coord(int i, int j, double* x, double* y) const = 0;

    void print_field() const;
};

class Pmesh : public Mesh {
public:
    Pmesh(int Nx, int Ny, double Lx, double Ly);
    double findx(int i) const override;
    double findy(int j) const override;
    void coord(int i, int j, double* x, double* y) const override;
};

class Umesh : public Mesh {
public:
    Umesh(int Nx, int Ny, double Lx, double Ly);
    double findx(int i) const override;
    double findy(int j) const override;
    void coord(int i, int j, double* x, double* y) const override;
};

class Vmesh : public Mesh {
public:
    Vmesh(int Nx, int Ny, double Lx, double Ly);
    double findx(int i) const override;
    double findy(int j) const override;
    void coord(int i, int j, double* x, double* y) const override;
};

#endif // MESH_H