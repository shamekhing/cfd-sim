#pragma once

#include "../core/Mesh.hpp"
#include "../core/Field.hpp"
#include "../core/UField.hpp"
#include "../core/VField.hpp"
#include "../core/Indexing.hpp"
#include "../linear/Vector.hpp"
#include "../linear/Matrix.hpp"

namespace assembly {

struct StaggeredMomentumSystem {
    linear::SparseMatrix A;
    linear::Vector b;
    linear::Vector faceDiagU;  // diagonal a_e per u-face (for d_e = A_e/a_e)
    linear::Vector faceDiagV;  // diagonal a_n per v-face (for d_n = A_n/a_n)
};

class StaggeredMomentumAssembler {
public:
    // Step I: Implicit momentum on faces. a_e u_e* = sum(a_nb u_nb*) + A_e(p_P* - p_E*) + S_u.
    // Steady: dt=0. Transient: pass dt and u_old; add (1/dt) to diagonal and u_old/dt to RHS.
    // source_u: optional body force in x-direction (e.g., dp/dx for channel flow)
    static StaggeredMomentumSystem assembleU(
        const core::Mesh& mesh,
        const core::UField& u,
        const core::VField& v,
        const core::Field& p,
        linear::Scalar nu,
        linear::Scalar dt = 0.0,
        const core::UField* u_old = nullptr,
        linear::Scalar source_u = 0.0)
    {
        using namespace core;
        using namespace linear;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int Nu = mesh.totalUFaces();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        linear::SparseBuilder builder(Nu, Nu);
        linear::Vector b(Nu);
        b.setZero();
        linear::Vector faceDiagU(Nu);
        faceDiagU.setZero();

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i <= Nx; ++i) {
                int row = idxU(i, j, Nx);
                if (isUBoundary(i, j, Nx, Ny)) {
                    builder.add(row, row, 1.0);
                    faceDiagU(row) = 1.0;
                    bool lid = (j == Ny - 1);
                    b(row) = lid ? 1.0 : 0.0;
                    continue;
                }
                // ── FVM coefficients (NOT per-unit-volume) ───────────────
                // All terms carry face-area or cell-volume scaling so that
                // the pressure gradient  b += A_e·(p_W - p_E)  is consistent.
                // This is essential when dx ≠ dy (non-square cells).
                linear::Scalar uC = u(i, j);
                linear::Scalar vFaceL = static_cast<linear::Scalar>(0.5) * (v(i - 1, j) + v(i - 1, j + 1));

                linear::Scalar aP = 0.0;
                linear::Scalar aE = 0.0, aW = 0.0, aN = 0.0, aS = 0.0;

                // Diffusion: ν·A_face / δ
                aE += nu * dy / dx;
                aW += nu * dy / dx;
                aN += nu * dx / dy;
                aS += nu * dx / dy;
                aP += 2.0 * nu * dy / dx + 2.0 * nu * dx / dy;

                // Upwind convection: F_face = u·A_face
                if (uC >= 0.0) {
                    aP += uC * dy;
                    aW += uC * dy;
                } else {
                    aP += (-uC) * dy;
                    aE += (-uC) * dy;
                }

                if (vFaceL >= 0.0) {
                    aP += vFaceL * dx;
                    aS += vFaceL * dx;
                } else {
                    aP += (-vFaceL) * dx;
                    aN += (-vFaceL) * dx;
                }

                // Transient: ΔV / Δt
                if (dt > 0.0) {
                    aP += dx * dy / dt;
                    b(row) += (u_old ? (*u_old)(i, j) : u(i, j)) * dx * dy / dt;
                }

                builder.add(row, row, aP);
                builder.add(row, idxU(i + 1, j, Nx), -aE);
                builder.add(row, idxU(i - 1, j, Nx), -aW);
                builder.add(row, idxU(i, j + 1, Nx), -aN);
                builder.add(row, idxU(i, j - 1, Nx), -aS);

                linear::Scalar pP = p(i - 1, j);
                linear::Scalar pE = p(i, j);
                linear::Scalar A_e = dy;
                b(row) += A_e * (pP - pE);
                
                // Add source term (body force) if provided
                if (source_u != 0.0) {
                    b(row) += source_u * dx * dy;  // Volume * body force
                }

                // ── Deferred correction: central differencing ────────────
                // The upwind scheme adds numerical diffusion ≈ |u|·Δx/2.
                // This explicit anti-diffusion source corrects toward a
                // central scheme while keeping the matrix diagonally dominant.
                // correction = –|F|·(φ_E + φ_W – 2φ_P)/2  per direction.
                {
                    linear::Scalar d2u_x = u(i + 1, j) - 2.0 * u(i, j) + u(i - 1, j);
                    linear::Scalar d2u_y = u(i, j + 1) - 2.0 * u(i, j) + u(i, j - 1);
                    b(row) -= std::abs(uC)     * dy * d2u_x * 0.5;
                    b(row) -= std::abs(vFaceL) * dx * d2u_y * 0.5;
                }

                faceDiagU(row) = aP;
            }
        }

        SparseMatrix A = builder.build();
        return {A, b, faceDiagU, linear::Vector(0)};
    }

    static StaggeredMomentumSystem assembleV(
        const core::Mesh& mesh,
        const core::UField& u,
        const core::VField& v,
        const core::Field& p,
        linear::Scalar nu,
        linear::Scalar dt = 0.0,
        const core::VField* v_old = nullptr,
        linear::Scalar source_v = 0.0)
    {
        using namespace core;
        using namespace linear;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int Nv = mesh.totalVFaces();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        SparseBuilder builder(Nv, Nv);
        Vector b(Nv);
        b.setZero();
        Vector faceDiagV(Nv);
        faceDiagV.setZero();

        for (int j = 0; j <= Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int row = idxV(i, j, Nx);
                if (isVBoundary(i, j, Nx, Ny)) {
                    builder.add(row, row, 1.0);
                    faceDiagV(row) = 1.0;
                    b(row) = 0.0;
                    continue;
                }
                // ── FVM coefficients (NOT per-unit-volume) ───────────────
                linear::Scalar vC = v(i, j);
                linear::Scalar uFaceB = static_cast<linear::Scalar>(0.5) * (u(i, j - 1) + u(i + 1, j - 1));

                linear::Scalar aP = 0.0;
                linear::Scalar aE = 0.0, aW = 0.0, aN = 0.0, aS = 0.0;

                // Diffusion: ν·A_face / δ
                aE += nu * dy / dx;
                aW += nu * dy / dx;
                aN += nu * dx / dy;
                aS += nu * dx / dy;
                aP += 2.0 * nu * dy / dx + 2.0 * nu * dx / dy;

                // Upwind convection: F_face = u·A_face
                if (uFaceB >= 0.0) {
                    aP += uFaceB * dy;
                    aW += uFaceB * dy;
                } else {
                    aP += (-uFaceB) * dy;
                    aE += (-uFaceB) * dy;
                }

                if (vC >= 0.0) {
                    aP += vC * dx;
                    aS += vC * dx;
                } else {
                    aP += (-vC) * dx;
                    aN += (-vC) * dx;
                }

                // Transient: ΔV / Δt
                if (dt > 0.0) {
                    aP += dx * dy / dt;
                    b(row) += (v_old ? (*v_old)(i, j) : v(i, j)) * dx * dy / dt;
                }

                builder.add(row, row, aP);
                builder.add(row, idxV(i + 1, j, Nx), -aE);
                builder.add(row, idxV(i - 1, j, Nx), -aW);
                builder.add(row, idxV(i, j + 1, Nx), -aN);
                builder.add(row, idxV(i, j - 1, Nx), -aS);

                linear::Scalar pP = p(i, j - 1);
                linear::Scalar pN = p(i, j);
                linear::Scalar A_n = dx;
                b(row) += A_n * (pP - pN);
                
                // Add source term (body force) if provided
                if (source_v != 0.0) {
                    b(row) += source_v * dx * dy;  // Volume * body force
                }

                // ── Deferred correction: central differencing ────────────
                {
                    linear::Scalar d2v_x = v(i + 1, j) - 2.0 * v(i, j) + v(i - 1, j);
                    linear::Scalar d2v_y = v(i, j + 1) - 2.0 * v(i, j) + v(i, j - 1);
                    b(row) -= std::abs(uFaceB) * dy * d2v_x * 0.5;
                    b(row) -= std::abs(vC)     * dx * d2v_y * 0.5;
                }

                faceDiagV(row) = aP;
            }
        }

        SparseMatrix A = builder.build();
        return {A, b, linear::Vector(0), faceDiagV};
    }
};

}
