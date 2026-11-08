
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
      m_up(Nx, Ny, Lx, Ly),
      m_vp(Nx, Ny, Lx, Ly),
      m_pseudoPressure(Nx, Ny, Lx, Ly)
{
    if (Nx < 2 || Ny < 2 || Lx <= 0.0 || Ly <= 0.0 || rho <= 0.0)
        throw std::invalid_argument("Invalid Solver dimensions or physical parameters");
    m_pseudoPressure.initialize([](double x, double y)
    {
	return 0;
    });


    m_up.initialize([](double x, double y)
    {
	return 0;
    });
    m_vp.initialize([](double x, double y)
    {
	return 0;
    });
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
    //Transforms the matrix index to flat index
    return static_cast<std::size_t>(j) * m_Nx + i;
}

void Solver::solvePoissonFFT()
{
    //Solve the poisson problem and fills the m_pseudoPressure vector with the solution 
    auto* rhs = static_cast<fftw_complex*>(m_rhs);
    auto* rhsHat = static_cast<fftw_complex*>(m_rhsHat);

    // Discrete divergence at pressure-cell centres. The mesh is periodic.
    for (int i = 0; i < m_Nx*m_Ny; ++i) {
        const double div = (m_up.getFieldValue(i) - m_up.getFieldValue(m_up.gowest(i))) / m_dx
                             + (m_vp.getFieldValue(i) - m_vp.getFieldValue(m_vp.gosouth(i))) / m_dy;
        rhs[i][0] = div;
        rhs[i][1] = 0.0;
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
    for (int i = 0; i < m_Nx*m_Ny; ++i)
    {
            m_pseudoPressure.setFieldValue(i, rhs[i][0] * normalization);
    }

}

// Override function solvePoissonFFT() to make the valoration of the case B with a given arbitrari u and v fields

void Solver::solvePoissonFFT(Umesh& u, Vmesh& v, Pmesh& p)
{
    auto* rhs = static_cast<fftw_complex*>(m_rhs);
    auto* rhsHat = static_cast<fftw_complex*>(m_rhsHat);

    // Discrete divergence at pressure-cell centres. The mesh is periodic.

    for (int i = 0; i< m_Nx*m_Ny; i++){
	const double div = (u.getFieldValue(i) - u.getFieldValue(u.gowest(i)))/m_dx + (v.getFieldValue(i) - v.getFieldValue(v.gosouth(i)))/m_dy;

	rhs[i][0] = div;
	rhs[i][1] = 0.0;
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
    for (int i = 0; i < m_Nx*m_Ny; ++i){
            m_pseudoPressure.setFieldValue(i, rhs[i][0] * normalization);
    }
    for (int i = 0; i < m_Nx*m_Ny; ++i) {
	u.setFieldValue(i, u.getFieldValue(i) - (m_pseudoPressure.getFieldValue(m_pseudoPressure.goeast(i)) - m_pseudoPressure.getFieldValue(i)) / m_dx);
    	v.setFieldValue(i, v.getFieldValue(i) - (m_pseudoPressure.getFieldValue(m_pseudoPressure.gonorth(i)) - m_pseudoPressure.getFieldValue(i)) / m_dy);
        p.setFieldValue(i,(m_rho / m_dt) * m_pseudoPressure.getFieldValue(i));
    }

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

    m_dtc = (velocity_max > 0.0)? m_dx / velocity_max : m_dtd;

    m_dt = 0.1 * std::min(m_dtd, m_dtc);

    // Do not advance beyond the requested final time.
    m_dt = std::min(m_dt, dtMax);

    // LOOP 1: convective/diffusive terms and AB2 predictor.
    // FVM functions are unchanged from the Part A verification.
   
    for (int i = 0; i < m_Nx*m_Ny; ++i) {
        const double convU = FVM::computeConvectiveU(i, u, v);
        const double difU  = FVM::computeDiffusiveU(i, u, m_nu);
        const double convV = FVM::computeConvectiveV(i, u, v);
        const double difV  = FVM::computeDiffusiveV(i, u, v, m_nu);

	// Momentum equation: du/dt = -C/V + D/V - (1/rho) grad(p).
        Ru[i] = (-convU + difU) / cellVolume;
        Rv[i] = (-convV + difV) / cellVolume;

	// Start AB2 with one forward-Euler step; using R^{-1}=0 in AB2
	// would incorrectly multiply the first increment by 3/2.
	if (m_step == 0) {
        	m_up.setFieldValue(i, u.getFieldValue(i) + m_dt * Ru[i]);
                m_vp.setFieldValue(i, v.getFieldValue(i) + m_dt * Rv[i]);
        } else {
                m_up.setFieldValue(i, u.getFieldValue(i) + m_dt * (1.5 * Ru[i] - 0.5 * m_RuPrevious[i]));
                m_vp.setFieldValue(i, v.getFieldValue(i) + m_dt * (1.5 * Rv[i] - 0.5 * m_RvPrevious[i]));
        }
    }


    solvePoissonFFT();

    // LOOP 2: correct velocities and recover physical pressure p = rho/dt * p_tilde.
    // u(i,j) lies between pressure cells i and i+1; v(i,j) between j and j+1.
    
    for (int i = 0; i < m_Nx*m_Ny; ++i) {
	u.setFieldValue(i, m_up.getFieldValue(i) - (m_pseudoPressure.getFieldValue(m_pseudoPressure.goeast(i)) - m_pseudoPressure.getFieldValue(i)) / m_dx);
    	v.setFieldValue(i, m_vp.getFieldValue(i) - (m_pseudoPressure.getFieldValue(m_pseudoPressure.gonorth(i)) - m_pseudoPressure.getFieldValue(i)) / m_dy);
        p.setFieldValue(i,(m_rho / m_dt) * m_pseudoPressure.getFieldValue(i));
    }

    // Keep R^n for the next AB2 step. Initial previous R is zero (R^{-1}=0).

    m_RuPrevious.swap(Ru);
    m_RvPrevious.swap(Rv);
    ++m_step;
    m_time += m_dt;

}

