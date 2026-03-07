#pragma once

#include <iostream>
#include <cmath>

#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "../include/assembly/MomentumAssembler.hpp"
#include "../include/linear/Solver.hpp"

namespace cases {

inline void runMomentumTest()
{
    const int Nx = 50;
    const int Ny = 50;

    core::Mesh mesh(Nx, Ny, 1.0, 1.0);

    core::Field u(mesh);
    core::Field v(mesh);
    core::Field p(mesh);   // dummy pressure for testing

    double nu = 0.01;

    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            u(i,j) = 1.0;
            v(i,j) = 0.0;
        }
    }

    const int maxIter = 20;

    for (int iter = 0; iter < maxIter; ++iter)
    {
        auto systemU =
            assembly::MomentumAssembler::assembleU(mesh, u, v, p, nu);

        linear::Vector newU =
            linear::BiCGSTABSolver().solve(systemU.A, systemU.b);
        u.raw() = newU;

        auto systemV =
            assembly::MomentumAssembler::assembleV(mesh, u, v, p, nu);

        linear::Vector newV =
            linear::BiCGSTABSolver().solve(systemV.A, systemV.b);
        v.raw() = newV;

        std::cout << "Iteration " << iter << " completed\n";
    }

    std::cout << "Momentum test finished\n";
}

}