#pragma once

#include "Field.hpp"

namespace core {

struct LidDrivenCavityBC {
    linear::Scalar u_lid = 1.0;
    int Nx = 0;
    int Ny = 0;
};

inline void applyLidDrivenCavityBC(Field& u, Field& v, const LidDrivenCavityBC& bc)
{
    const int Nx = bc.Nx;
    const int Ny = bc.Ny;

    for (int j = 0; j < Ny; ++j) {
        u(0, j) = 0.0;
        u(Nx - 1, j) = 0.0;
    }
    for (int i = 0; i < Nx; ++i) {
        u(i, 0) = 0.0;
        u(i, Ny - 1) = bc.u_lid;
    }

    for (int j = 0; j < Ny; ++j) {
        v(0, j) = 0.0;
        v(Nx - 1, j) = 0.0;
    }
    for (int i = 0; i < Nx; ++i) {
        v(i, 0) = 0.0;
        v(i, Ny - 1) = 0.0;
    }
}

}
