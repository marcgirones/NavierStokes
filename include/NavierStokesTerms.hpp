#ifndef NAVIER_STOKES_TERMS_H
#define NAVIER_STOKES_TERMS_H

#include "Mesh.hpp"

namespace FVM {

    // Términos para el componente U de la velocidad
    double computeConvectiveU(int i, const Umesh& u, const Vmesh& v);
    double computeDiffusiveU(int i, const Umesh& u, double nu);


    // Términos para el componente V de la velocidad
    double computeConvectiveV(int i, const Umesh& u, const Vmesh& v);
    double computeDiffusiveV(int i, const Umesh& u, const Vmesh& v, double nu);

} // namespace FVM

#endif // NAVIER_STOKES_TERMS_H
