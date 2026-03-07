#pragma once

#include "../core/Mesh.hpp"
#include "../core/Field.hpp"
#include "../core/Indexing.hpp"
#include "../linear/Vector.hpp"
#include "../linear/Matrix.hpp"

namespace assembly {

struct PressureSystem {
    linear::SparseMatrix A;
    linear::Vector b;
};

class PressureCorrectionAssembler {
public:

    static PressureSystem assemble(
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v)
    {
        using namespace core;
        using namespace linear;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int N  = mesh.totalCells();

        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        SparseBuilder builder(N, N);
        Vector b(N);
        b.setZero();

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                int P = id(i, j, Nx);

                if (mesh.isBoundary(i,j))
                {
                    builder.add(P, P, 1.0);
                    b(P) = 0.0;
                    continue;
                }

                int E = id(i+1, j, Nx);
                int W = id(i-1, j, Nx);
                int Nn = id(i, j+1, Nx);
                int S = id(i, j-1, Nx);

                linear::Scalar aP = -2.0/(dx*dx) - 2.0/(dy*dy);
                linear::Scalar aE =  1.0/(dx*dx);
                linear::Scalar aW =  1.0/(dx*dx);
                linear::Scalar aN =  1.0/(dy*dy);
                linear::Scalar aS =  1.0/(dy*dy);

                builder.add(P, P, aP);
                builder.add(P, E, aE);
                builder.add(P, W, aW);
                builder.add(P, Nn, aN);
                builder.add(P, S, aS);

                                linear::Scalar div =
                    (u(i+1,j) - u(i-1,j))/(2.0*dx)
                  + (v(i,j+1) - v(i,j-1))/(2.0*dy);

                b(P) = div;
            }
        }

        SparseMatrix A = builder.build();
        return {A, b};
    }
};

}