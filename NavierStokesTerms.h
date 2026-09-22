#ifndef NAVIER_STOKES_TERMS_H
#define NAVIER_STOKES_TERMS_H

#include "Mesh.h"

namespace FVM {

    // Términos para el componente U de la velocidad
    double computeConvectiveU(int i, int j, const Umesh& u, const Vmesh& v);
    double computeDiffusiveU(int i, int j, const Umesh& u, double nu);

    // Términos para el componente V de la velocidad
    double computeConvectiveV(int i, int j, const Umesh& u, const Vmesh& v);
    double computeDiffusiveV(int i, int j, const Umesh& u, const Vmesh& v, double nu);

} // namespace FVM

#endif // NAVIER_STOKES_TERMS_H