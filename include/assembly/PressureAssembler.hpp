#pragma once

#include "../core/Mesh.hpp"
#include "../core/Field.hpp"
#include "../core/Indexing.hpp"
#include "../linear/Vector.hpp"
#include "../linear/Matrix.hpp"

namespace assembly {

struct PoissonSystem {
    linear::SparseMatrix A;
    linear::Vector b;
};

class PressureAssembler {
public:

    static PoissonSystem assemble(
        const core::Mesh& mesh,
        const core::Field& rhs,
        linear::Scalar dirichletValue = 0.0)
    {
        using namespace linear;
        using namespace core;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int N  = mesh.totalCells();

        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        const linear::Scalar dx2 = dx * dx;
        const linear::Scalar dy2 = dy * dy;

        SparseBuilder builder(N, N);
        Vector b(N);
        b.setZero();

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                int P = id(i, j, Nx);

                if (mesh.isBoundary(i, j))
                {
                    // Dirichlet boundary condition
                    builder.add(P, P, 1.0);
                    b(P) = dirichletValue;
                }
                else
                {
                    int E = id(i+1, j, Nx);
                    int W = id(i-1, j, Nx);
                    int Nn = id(i, j+1, Nx);
                    int S = id(i, j-1, Nx);

                    linear::Scalar aP = -2.0/dx2 - 2.0/dy2;
                    linear::Scalar aE =  1.0/dx2;
                    linear::Scalar aW =  1.0/dx2;
                    linear::Scalar aN =  1.0/dy2;
                    linear::Scalar aS =  1.0/dy2;

                    builder.add(P, P, aP);
                    builder.add(P, E, aE);
                    builder.add(P, W, aW);
                    builder.add(P, Nn, aN);
                    builder.add(P, S, aS);

                    b(P) = rhs(i, j);
                }
            }
        }

        SparseMatrix A = builder.build();

        return {A, b};
    }
};

}