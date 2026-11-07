#include "NavierStokesTerms.hpp"

namespace FVM {

// ============================================================================
// U Convection
// Integrated term: ∫ (∇ · (u u)) dV ≈ (F_e*u_e - F_w*u_w) + (F_n*u_n - F_s*u_s)
// ============================================================================
double computeConvectiveU(int i, const Umesh& u, const Vmesh& v) {
    const double dy = u.getLey();
    const double dx = u.getLex();

    // 1. Values transported through the faces (Average centered in U)
    double u_e = 0.5 * (u.getFieldValue(u.goeast(i)) + u.getFieldValue(i));
    double u_w = 0.5 * (u.getFieldValue(u.gowest(i)) + u.getFieldValue(i));
    double u_n = 0.5 * (u.getFieldValue(u.gonorth(i)) + u.getFieldValue(i));
    double u_s = 0.5 * (u.getFieldValue(u.gosouth(i)) + u.getFieldValue(i));

    // 2. Mass flow rate computed in its natural position
    // Faces East/West (U-velocity)
    double F_e = 0.5 * (u.getFieldValue(u.goeast(i)) + u.getFieldValue(i)) * dy;
    double F_w = 0.5 * (u.getFieldValue(u.gowest(i)) + u.getFieldValue(i)) * dy;

    // Faces North/South (V-velocity interpolated in the U volume)
    double F_n = 0.5 * (v.getFieldValue(i) + v.getFieldValue(v.goeast(i))) * dx;
    double F_s = 0.5 * (v.getFieldValue(v.gosouth(i)) + v.getFieldValue(v.gosouth(v.goeast(i)))) * dx;

    return (F_e * u_e - F_w * u_w) + (F_n * u_n - F_s * u_s);
}
// ============================================================================
// DIFUSIÓN PARA U
// U Diffusion
// Integrated term: ∫ (∇ · (ν ∇u)) dV ≈ ν * [ Δy*(∂u/∂x)_e - Δy*(∂u/∂x)_w 
//                                           + Δx*(∂u/∂y)_n - Δx*(∂u/∂y)_s ]
// ============================================================================

double computeDiffusiveU(int i, const Umesh& u, double nu) {
    const double dx = u.getLex();
    const double dy = u.getLey();

    // Face derivatives of the control volume U
    double dudx_e = (u.getFieldValue(u.goeast(i)) - u.getFieldValue(i)) / dx;
    double dudx_w = (u.getFieldValue(i) - u.getFieldValue(u.gowest(i))) / dx;
    double dudy_n = (u.getFieldValue(u.gonorth(i)) - u.getFieldValue(i))     / dy;
    double dudy_s = (u.getFieldValue(i) - u.getFieldValue(u.gosouth(i))) / dy;

    return nu * ((dudx_e - dudx_w) * dy + (dudy_n - dudy_s) * dx);
}
// ============================================================================
// V Conveccion
// ============================================================================

double computeConvectiveV(int i, const Umesh& u, const Vmesh& v) {
    const double dx = v.getLex();
    const double dy = v.getLey();

    // Transported values of V in the faces
    double v_e = 0.5 * (v.getFieldValue(v.goeast(i)) + v.getFieldValue(i));
    double v_w = 0.5 * (v.getFieldValue(v.gowest(i)) + v.getFieldValue(i));
    double v_n = 0.5 * (v.getFieldValue(v.gonorth(i)) + v.getFieldValue(i));
    double v_s = 0.5 * (v.getFieldValue(v.gosouth(i)) + v.getFieldValue(i));

    // Mass flow rate interpolated in the V volume
    double F_e = 0.5 * (u.getFieldValue(i) + u.getFieldValue(u.gonorth(i))) * dy;
    double F_w = 0.5 * (u.getFieldValue(u.gowest(i)) + u.getFieldValue(u.gonorth(u.gowest(i)))) * dy;
    double F_n = 0.5 * (v.getFieldValue(v.gonorth(i)) + v.getFieldValue(i)) * dx;
    double F_s = 0.5 * (v.getFieldValue(i) + v.getFieldValue(v.gosouth(i))) * dx;

    return (F_e * v_e - F_w * v_w) + (F_n * v_n - F_s * v_s);
}
// ============================================================================
// V Diffusion
// ============================================================================

double computeDiffusiveV(int i, const Umesh& u, const Vmesh& v, double nu) {
    const double dx = v.getLex();
    const double dy = v.getLey();

    double dvdx_e = (v.getFieldValue(v.goeast(i)) - v.getFieldValue(i)) / dx;
    double dvdx_w = (v.getFieldValue(i) - v.getFieldValue(v.gowest(i))) / dx;
    double dvdy_n = (v.getFieldValue(v.gonorth(i)) - v.getFieldValue(i)) / dy;
    double dvdy_s = (v.getFieldValue(i) - v.getFieldValue(v.gosouth(i))) / dy;

    return nu * ((dvdx_e - dvdx_w) * dy + (dvdy_n - dvdy_s) * dx);
}

} // namespace FVM
