#pragma once

#include <iostream>

#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "../include/assembly/ScalarAssembler.hpp"
#include "../include/linear/Solver.hpp"

namespace cases {

inline void runScalarCase()
{
    const int Nx = 80;
    const int Ny = 40;

    core::Mesh mesh(Nx, Ny, 1.0, 1.0);

    double u = 1.0;
    double gamma = 0.01;

    auto system =
        assembly::ScalarAssembler::assemble(
            mesh,
            u,
            gamma,
            assembly::ConvectionScheme::Upwind);

    linear::Vector phi =
        linear::BiCGSTABSolver().solve(
            system.A,
            system.b);

    std::cout << "Scalar Convection-Diffusion Case Completed\n";
}

}  // namespace cases

