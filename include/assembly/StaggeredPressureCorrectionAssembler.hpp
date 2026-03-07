#pragma once

#include "../core/Mesh.hpp"
#include "../core/Field.hpp"
#include "../core/UField.hpp"
#include "../core/VField.hpp"
#include "../core/Indexing.hpp"
#include "../linear/Vector.hpp"
#include "../linear/Matrix.hpp"

namespace assembly {

struct StaggeredPressureSystem {
    linear::SparseMatrix A;
    linear::Vector b;
};

// Step III: a_P p_P' = a_E p_E' + a_W p_W' + a_S p_S' + a_N p_N' + b
// a_E = d_e A_e, d_e = A_e/a_e from momentum; b = (u_w* - u_e*)A + (v_s* - v_n*)A
class StaggeredPressureCorrectionAssembler {
public:
    static StaggeredPressureSystem assemble(
        const core::Mesh& mesh,
        const core::UField& ustar,
        const core::VField& vstar,
        const linear::Vector& faceDiagU,
        const linear::Vector& faceDiagV)
    {
        using namespace core;
        using namespace linear;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int N = mesh.totalCells();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        SparseBuilder builder(N, N);
        Vector b(N);
        b.setZero();

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int P = id(i, j, Nx);
                if (mesh.isBoundary(i, j)) {
                    // Well-posed pressure-correction BC for this implementation:
                    // p' = 0 on all outer boundaries.
                    builder.add(P, P, 1.0);
                    b(P) = 0.0;
                    continue;
                }
                int E = id(i + 1, j, Nx);
                int W = id(i - 1, j, Nx);
                int Nn = id(i, j + 1, Nx);
                int S = id(i, j - 1, Nx);

                linear::Scalar a_e = faceDiagU(idxU(i + 1, j, Nx));
                linear::Scalar a_w = faceDiagU(idxU(i, j, Nx));
                linear::Scalar a_n = faceDiagV(idxV(i, j + 1, Nx));
                linear::Scalar a_s = faceDiagV(idxV(i, j, Nx));

                linear::Scalar A_e = dy, A_w = dy, A_n = dx, A_s = dx;
                linear::Scalar d_e = (a_e > static_cast<linear::Scalar>(1e-30)) ? A_e / a_e : static_cast<linear::Scalar>(0.0);
                linear::Scalar d_w = (a_w > static_cast<linear::Scalar>(1e-30)) ? A_w / a_w : static_cast<linear::Scalar>(0.0);
                linear::Scalar d_n = (a_n > static_cast<linear::Scalar>(1e-30)) ? A_n / a_n : static_cast<linear::Scalar>(0.0);
                linear::Scalar d_s = (a_s > static_cast<linear::Scalar>(1e-30)) ? A_s / a_s : static_cast<linear::Scalar>(0.0);

                linear::Scalar aE = d_e * A_e;
                linear::Scalar aW = d_w * A_w;
                linear::Scalar aN = d_n * A_n;
                linear::Scalar aS = d_s * A_s;
                linear::Scalar aP = aE + aW + aN + aS;

                // Positive diagonal, negative off-diagonals (Poisson-like form)
                builder.add(P, P, aP);
                builder.add(P, E, -aE);
                builder.add(P, W, -aW);
                builder.add(P, Nn, -aN);
                builder.add(P, S, -aS);

                linear::Scalar u_w = ustar(i, j);
                linear::Scalar u_e = ustar(i + 1, j);
                linear::Scalar v_s = vstar(i, j);
                linear::Scalar v_n = vstar(i, j + 1);
                b(P) = (u_w - u_e) * dy + (v_s - v_n) * dx;
            }
        }

        SparseMatrix A = builder.build();
        return {A, b};
    }
};

}
