#pragma once

#include "../core/Mesh.hpp"
#include "../core/Field.hpp"
#include "../core/Indexing.hpp"
#include "../linear/Vector.hpp"
#include "../linear/Matrix.hpp"

namespace assembly {

enum class ConvectionScheme {
    Central,
    Upwind
};

struct ScalarSystem {
    linear::SparseMatrix A;
    linear::Vector b;
};

class ScalarAssembler {
public:

    static ScalarSystem assemble(
        const core::Mesh& mesh,
        linear::Scalar u,
        linear::Scalar gamma,
        ConvectionScheme scheme)
    {
        using namespace core;
        using namespace linear;

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

                if (mesh.isLeftBoundary(i))
                {
                    builder.add(P, P, 1.0);
                    b(P) = 1.0;
                    continue;
                }

                if (mesh.isRightBoundary(i) ||
                    mesh.isTopBoundary(j) ||
                    mesh.isBottomBoundary(j))
                {
                    builder.add(P, P, 1.0);
                    b(P) = 0.0;
                    continue;
                }

                int E = id(i+1, j, Nx);
                int W = id(i-1, j, Nx);
                int Nn = id(i, j+1, Nx);
                int S = id(i, j-1, Nx);

                linear::Scalar aP = 0.0;
                linear::Scalar aE = 0.0;
                linear::Scalar aW = 0.0;
                linear::Scalar aN = 0.0;
                linear::Scalar aS = 0.0;

                // Diffusion
                aE += gamma/dx2;
                aW += gamma/dx2;
                aN += gamma/dy2;
                aS += gamma/dy2;
                aP -= 2.0*gamma/dx2 + 2.0*gamma/dy2;

                // Convection
                if (scheme == ConvectionScheme::Central)
                {
                    aE -= u/(2.0*dx);
                    aW += u/(2.0*dx);
                }
                else if (scheme == ConvectionScheme::Upwind)
                {
                    if (u > 0)
                    {
                        aP += u/dx;
                        aW -= u/dx;
                    }
                    else
                    {
                        aP -= u/dx;
                        aE += u/dx;
                    }
                }

                builder.add(P, P, aP);
                builder.add(P, E, aE);
                builder.add(P, W, aW);
                builder.add(P, Nn, aN);
                builder.add(P, S, aS);
            }
        }

        SparseMatrix A = builder.build();
        return {A, b};
    }
};

}