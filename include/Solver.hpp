
#ifndef SOLVER_H
#define SOLVER_H

#include "Mesh.hpp"
#include <vector>

// Time integrator and pressure-velocity coupling for the periodic 2-D case.
// The convective/diffusive discretizations remain in namespace FVM.
class Solver {
public:
    Solver(int Nx, int Ny, double Lx, double Ly,
           double nu, double rho);
    ~Solver();

    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;

    // One complete time step: AB2 predictor -> Poisson solve -> correction.
    void step(Umesh& u, Vmesh& v, Pmesh& p, double dtMax);

    void solvePoissonFFT(Umesh& u, Vmesh& v, Pmesh& p);

    int getStep() const { return m_step; }
    double getTime() const { return m_time; }
    double getdt() const { return m_dt; }

private:
    int m_Nx, m_Ny;
    double m_Lx, m_Ly, m_dx, m_dy;
    double m_nu, m_rho, m_dtd, m_dtc, m_dt;
    int m_step = 0;
    double m_time = 0.0;

    // Previous R^n fields, flattened as j*Nx+i. R^{-1}=0 initially.
    std::vector<double> m_RuPrevious, m_RvPrevious;

    // Predictor velocities and pseudopressure are stored separately so that
    // updating one cell never changes data needed by another cell.
    std::vector<double> m_up, m_vp, m_pseudoPressure;

    // FFTW resources are hidden from the header using void pointers.
    void* m_rhs = nullptr;
    void* m_rhsHat = nullptr;
    void* m_forwardPlan = nullptr;
    void* m_inversePlan = nullptr;

    std::size_t index(int i, int j) const;
    void solvePoissonFFT();
};

#endif
