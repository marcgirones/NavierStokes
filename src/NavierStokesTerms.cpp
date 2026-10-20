#include "NavierStokesTerms.hpp"

namespace FVM {

// ============================================================================
// CONVECCIÓN PARA U
// Integrated term: ∫ (∇ · (u u)) dV ≈ (F_e*u_e - F_w*u_w) + (F_n*u_n - F_s*u_s)
// ============================================================================
double computeConvectiveU(int i, int j, const Umesh& u, const Vmesh& v) {
    const double dy = u.getLey();
    const double dx = u.getLex();

    // 1. Valores transportados en las caras (Promedio centrado de U)
    double u_e = 0.5 * (u(i + 1, j) + u(i, j));
    double u_w = 0.5 * (u(i - 1, j) + u(i, j));
    double u_n = 0.5 * (u(i, j + 1) + u(i, j));
    double u_s = 0.5 * (u(i, j - 1) + u(i, j));

    // 2. Flujos masicos evaluados en sus posiciones naturales
    // Caras Este/Oeste (U-velocidad)
    double F_e = 0.5 * (u(i + 1, j) + u(i, j)) * dy;
    double F_w = 0.5 * (u(i - 1, j) + u(i, j)) * dy;

    // Caras Norte/Sur (V-velocidad interpolada al volumen de U)
    double F_n = 0.5 * (v(i, j)     + v(i + 1, j))     * dx;
    double F_s = 0.5 * (v(i, j - 1) + v(i + 1, j - 1)) * dx;

    return (F_e * u_e - F_w * u_w) + (F_n * u_n - F_s * u_s);
}

// ============================================================================
// DIFUSIÓN PARA U
// Integrated term: ∫ (∇ · (ν ∇u)) dV ≈ ν * [ Δy*(∂u/∂x)_e - Δy*(∂u/∂x)_w 
//                                           + Δx*(∂u/∂y)_n - Δx*(∂u/∂y)_s ]
// ============================================================================
double computeDiffusiveU(int i, int j, const Umesh& u, double nu) {
    const double dx = u.getLex();
    const double dy = u.getLey();

    // Derivadas en las caras del volumen de control de U
    double dudx_e = (u(i + 1, j) - u(i, j))     / dx;
    double dudx_w = (u(i, j)     - u(i - 1, j)) / dx;
    double dudy_n = (u(i, j + 1) - u(i, j))     / dy;
    double dudy_s = (u(i, j)     - u(i, j - 1)) / dy;

    return nu * ((dudx_e - dudx_w) * dy + (dudy_n - dudy_s) * dx);
}

// ============================================================================
// CONVECCIÓN PARA V
// ============================================================================
double computeConvectiveV(int i, int j, const Umesh& u, const Vmesh& v) {
    const double dx = v.getLex();
    const double dy = v.getLey();

    // Valores transportados V en las caras
    double v_e = 0.5 * (v(i + 1, j) + v(i, j));
    double v_w = 0.5 * (v(i - 1, j) + v(i, j));
    double v_n = 0.5 * (v(i, j + 1) + v(i, j));
    double v_s = 0.5 * (v(i, j - 1) + v(i, j));

    // Flujos masicos interpolados al volumen de V
    double F_e = 0.5 * (u(i, j)     + u(i, j + 1))     * dy;
    double F_w = 0.5 * (u(i - 1, j) + u(i - 1, j + 1)) * dy;
    double F_n = 0.5 * (v(i, j + 1) + v(i, j))         * dx;
    double F_s = 0.5 * (v(i, j - 1) + v(i, j))         * dx;

    return (F_e * v_e - F_w * v_w) + (F_n * v_n - F_s * v_s);
}

// ============================================================================
// DIFUSIÓN PARA V
// ============================================================================
double computeDiffusiveV(int i, int j, const Umesh& u, const Vmesh& v, double nu) {
    const double dx = v.getLex();
    const double dy = v.getLey();

    double dvdx_e = (v(i + 1, j) - v(i, j))     / dx;
    double dvdx_w = (v(i, j)     - v(i - 1, j)) / dx;
    double dvdy_n = (v(i, j + 1) - v(i, j))     / dy;
    double dvdy_s = (v(i, j)     - v(i, j - 1)) / dy;

    return nu * ((dvdx_e - dvdx_w) * dy + (dvdy_n - dvdy_s) * dx);
}

} // namespace FVM
