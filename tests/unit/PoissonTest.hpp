#pragma once

#include <iostream>
#include <cmath>

#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "../include/core/Indexing.hpp"
#include "../include/assembly/PressureAssembler.hpp"
#include "../include/linear/Solver.hpp"

namespace cases {

inline void runPoissonTest()
{
    const int Nx = 50;
    const int Ny = 50;

    core::Mesh mesh(Nx, Ny, 1.0, 1.0);
    core::Field rhs(mesh);

    // Manufactured solution:
    // phi = sin(pi x) sin(pi y)
    // Laplacian(phi) = -2 pi^2 sin(pi x) sin(pi y)

    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            double x = mesh.x(i);
            double y = mesh.y(j);

            rhs(i,j) =
                -2.0 * M_PI * M_PI *
                std::sin(M_PI * x) *
                std::sin(M_PI * y);
        }
    }

    auto system =
        assembly::PressureAssembler::assemble(mesh, rhs, 0.0);

    linear::Vector phi =
        linear::BiCGSTABSolver().solve(system.A, system.b);

    double error = 0.0;

    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            double x = mesh.x(i);
            double y = mesh.y(j);

            double exact =
                std::sin(M_PI * x) *
                std::sin(M_PI * y);

            int P = core::id(i, j, Nx);

            error += std::pow(phi(P) - exact, 2);
        }
    }

    error = std::sqrt(error / (Nx * Ny));

    std::cout << "Poisson Test Completed\n";
    std::cout << "Grid: " << Nx << " x " << Ny << "\n";
    std::cout << "L2 Error: " << error << "\n";
}

}  // namespace cases

