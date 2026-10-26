
#include "Solver.hpp"
#include "NavierStokesTerms.hpp"

#include <fftw3.h>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <chrono>

Solver::Solver(int Nx, int Ny, double Lx, double Ly,
               double nu, double rho)
    : m_Nx(Nx), m_Ny(Ny), m_Lx(Lx), m_Ly(Ly),
      m_dx(Lx / Nx), m_dy(Ly / Ny), m_nu(nu), m_rho(rho),
      m_dtd(0.5*((Lx/Nx)*(Lx/Nx))/nu), m_dtc(m_dtd), m_dt(m_dtc),
      m_RuPrevious(static_cast<std::size_t>(Nx) * Ny, 0.0),
      m_RvPrevious(static_cast<std::size_t>(Nx) * Ny, 0.0),
      m_up(static_cast<std::size_t>(Nx) * Ny, 0.0),
      m_vp(static_cast<std::size_t>(Nx) * Ny, 0.0),
      m_pseudoPressure(static_cast<std::size_t>(Nx) * Ny, 0.0)
{
    if (Nx < 2 || Ny < 2 || Lx <= 0.0 || Ly <= 0.0 || rho <= 0.0)
        throw std::invalid_argument("Invalid Solver dimensions or physical parameters");

    const std::size_t n = static_cast<std::size_t>(Nx) * Ny;
    auto* rhs = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));
    auto* rhsHat = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));
    if (!rhs || !rhsHat) {
        if (rhs) fftw_free(rhs);
        if (rhsHat) fftw_free(rhsHat);
        throw std::bad_alloc();
    }
    m_rhs = rhs;
    m_rhsHat = rhsHat;
    m_forwardPlan = fftw_plan_dft_2d(m_Ny, m_Nx, rhs, rhsHat,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
    m_inversePlan = fftw_plan_dft_2d(m_Ny, m_Nx, rhsHat, rhs,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
    if (!m_forwardPlan || !m_inversePlan)
        throw std::runtime_error("FFTW could not create the 2-D transform plans");
}

Solver::~Solver()
{
    if (m_forwardPlan) fftw_destroy_plan(static_cast<fftw_plan>(m_forwardPlan));
    if (m_inversePlan) fftw_destroy_plan(static_cast<fftw_plan>(m_inversePlan));
    if (m_rhs) fftw_free(m_rhs);
    if (m_rhsHat) fftw_free(m_rhsHat);
}

std::size_t Solver::index(int i, int j) const
{
    return static_cast<std::size_t>(j) * m_Nx + i;
}

void Solver::solvePoissonFFT(const Umesh& u, const Vmesh& v)
{
    auto* rhs = static_cast<fftw_complex*>(m_rhs);
    auto* rhsHat = static_cast<fftw_complex*>(m_rhsHat);

    // Discrete divergence at pressure-cell centres. The mesh is periodic.
    for (int j = 0; j < m_Ny; ++j) {
        for (int i = 0; i < m_Nx; ++i) {
            const double div = (m_up[index(i, j)] - m_up[index(i - 1 < 0 ? m_Nx - 1 : i - 1, j)]) / m_dx
                             + (m_vp[index(i, j)] - m_vp[index(i, j - 1 < 0 ? m_Ny - 1 : j - 1)]) / m_dy;
            rhs[index(i, j)][0] = div;
            rhs[index(i, j)][1] = 0.0;
        }
    }

    fftw_execute(static_cast<fftw_plan>(m_forwardPlan));

    // Invert the eigenvalues of the same second-order discrete Laplacian.
    // The (kx,ky)=(0,0) mode is singular: set mean pseudopressure to zero.
    for (int ky = 0; ky < m_Ny; ++ky) {
        const double sy = std::sin(M_PI * ky / m_Ny);
        for (int kx = 0; kx < m_Nx; ++kx) {
            const std::size_t k = index(kx, ky);
            if (kx == 0 && ky == 0) {
                rhsHat[k][0] = 0.0;
                rhsHat[k][1] = 0.0;
                continue;
            }
            const double sx = std::sin(M_PI * kx / m_Nx);
            const double lambda = -4.0 * sx * sx / (m_dx * m_dx)
                                  -4.0 * sy * sy / (m_dy * m_dy);
            rhsHat[k][0] /= lambda;
            rhsHat[k][1] /= lambda;
        }
    }

    fftw_execute(static_cast<fftw_plan>(m_inversePlan));
    const double normalization = 1.0 / (static_cast<double>(m_Nx) * m_Ny);
    for (int j = 0; j < m_Ny; ++j)
        for (int i = 0; i < m_Nx; ++i)
            m_pseudoPressure[index(i, j)] = rhs[index(i, j)][0] * normalization;

    (void)u;
    (void)v;
}

void Solver::step(Umesh& u, Vmesh& v, Pmesh& p, double dtMax)
{
    //auto start_step = std::chrono::high_resolution_clock::now();
    const double cellVolume = m_dx * m_dy;
    std::vector<double> Ru(m_RuPrevious.size());
    std::vector<double> Rv(m_RvPrevious.size());

    // Determine the timestep once, before computing the predictor.
    const double u_max = u.max_val();
    const double v_max = v.max_val();
    const double velocity_max = std::max(u_max, v_max);

    m_dtc = (velocity_max > 0.0)
          ? m_dx / velocity_max
          : m_dtd;

    m_dt = 0.1 * std::min(m_dtd, m_dtc);

    // Do not advance beyond the requested final time.
    m_dt = std::min(m_dt, dtMax);

    // LOOP 1: convective/diffusive terms and AB2 predictor.
    // FVM functions are unchanged from the Part A verification.
    
    //auto start_cd = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < m_Nx; ++i) {
        for (int j = 0; j < m_Ny; ++j) {
            const std::size_t k = index(i, j);
            const double convU = FVM::computeConvectiveU(i, j, u, v);
            const double difU  = FVM::computeDiffusiveU(i, j, u, m_nu);
            const double convV = FVM::computeConvectiveV(i, j, u, v);
            const double difV  = FVM::computeDiffusiveV(i, j, u, v, m_nu);

            // Momentum equation: du/dt = -C/V + D/V - (1/rho) grad(p).
            Ru[k] = (-convU + difU) / cellVolume;
            Rv[k] = (-convV + difV) / cellVolume;

            // Start AB2 with one forward-Euler step; using R^{-1}=0 in AB2
            // would incorrectly multiply the first increment by 3/2.
            if (m_step == 0) {
                m_up[k] = u(i, j) + m_dt * Ru[k];
                m_vp[k] = v(i, j) + m_dt * Rv[k];
            } else {
                m_up[k] = u(i, j) + m_dt * (1.5 * Ru[k] - 0.5 * m_RuPrevious[k]);
                m_vp[k] = v(i, j) + m_dt * (1.5 * Rv[k] - 0.5 * m_RvPrevious[k]);
            }
        }
    }

    //auto end_cd = std::chrono::high_resolution_clock::now();
    // Solve Laplacian(pseudoPressure) = div(predictor velocity).
    //auto start_p = std::chrono::high_resolution_clock::now();

    solvePoissonFFT(u, v);

    //auto end_p = std::chrono::high_resolution_clock::now();
    // LOOP 2: correct velocities and recover physical pressure p = rho/dt * p_tilde.
    // u(i,j) lies between pressure cells i and i+1; v(i,j) between j and j+1.
    
    double maxgrad =0;
    //auto start_correction = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < m_Nx; ++i) {
        for (int j = 0; j < m_Ny; ++j) {
            const std::size_t k = index(i, j);
            const int ip = (i + 1) % m_Nx;
            const int jp = (j + 1) % m_Ny;
            u(i, j) = m_up[k] - (m_pseudoPressure[index(ip, j)] - m_pseudoPressure[k]) / m_dx;
            v(i, j) = m_vp[k] - (m_pseudoPressure[index(i, jp)] - m_pseudoPressure[k]) / m_dy;
            p(i, j) = (m_rho / m_dt) * m_pseudoPressure[k];
	    double gradu = u.gradi(i);
	    double gradv = v.gradi(i);
	    double grad = std::abs(gradu+gradv);
	    if(grad>maxgrad) maxgrad = grad;
   
	}
    }


    /*for (int i = 0; i < m_Nx*m_Ny; ++i) {
	u.setFieldValue(i, m_up[i] - (m_pseudoPressure[u.goeast(i)] - m_pseudoPressure[i]) / m_dx);
    	v.setFieldValue(i, m_vp[i] - (m_pseudoPressure[v.gonorth(i)] - m_pseudoPressure[i]) / m_dy);
        p.setFieldValue(i,(m_rho / m_dt) * m_pseudoPressure[i]);
	double gradu = u.gradi(i);
	double gradv = v.gradi(i);
	double grad = std::abs(gradu+gradv);
	if(grad>maxgrad) maxgrad = grad;
    }*/

    std::cout<<"gradient màxim de u^(n+1) "<<maxgrad<<"\n";
    //auto end_correction = std::chrono::high_resolution_clock::now();
    // Keep R^n for the next AB2 step. Initial previous R is zero (R^{-1}=0).

    //auto start_update = std::chrono::high_resolution_clock::now();
    m_RuPrevious.swap(Ru);
    m_RvPrevious.swap(Rv);
    ++m_step;
    m_time += m_dt;

    //auto end_update = std::chrono::high_resolution_clock::now();
    //auto end_step = std::chrono::high_resolution_clock::now();
    

    //auto duration_step = std::chrono::duration_cast<std::chrono::microseconds>(end_step - start_step);
    //auto duration_cd = std::chrono::duration_cast<std::chrono::microseconds>(end_cd - start_cd);
    //auto duration_p = std::chrono::duration_cast<std::chrono::microseconds>(end_p - start_p);
    //auto duration_correction = std::chrono::duration_cast<std::chrono::microseconds>(end_correction - start_correction);
    //auto duration_update = std::chrono::duration_cast<std::chrono::microseconds>(end_update - start_update);

    //std::cout<<"Duration: Total episode --> "<<duration_step.count()<<" || Conv/dif part --> "<<duration_cd.count()<<" || Pressure part --> "<<duration_p.count()<<" || Correction --> "<<duration_correction.count()<<" || Update --> "<<duration_update.count()<<"\n";
}

