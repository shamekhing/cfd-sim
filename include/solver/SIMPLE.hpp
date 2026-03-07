#pragma once

#include "../include/linear/BiCGSTAB.hpp"
#include "../include/linear/LU.hpp"
#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "../include/assembly/MomentumAssembler.hpp"
#include "../include/assembly/PressureCorrectionAssembler.hpp"

namespace solver {

class SIMPLE {
public:
    // Explicit momentum predictor (reference-style: u* = u_old + dt*(-conv + nu*laplacian - grad p))
    // Uses current u,v (with BC) for convection/diffusion; lid influence enters via stencil.
    static void explicitMomentumPredictor(
        const core::Mesh& mesh,
        const core::Field& u_old,
        const core::Field& v_old,
        const core::Field& u,
        const core::Field& v,
        const core::Field& p,
        linear::Scalar nu,
        linear::Scalar dt,
        core::Field& ustar,
        core::Field& vstar)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();
        const linear::Scalar dx2 = dx * dx;
        const linear::Scalar dy2 = dy * dy;

        ustar.raw() = u_old.raw();
        vstar.raw() = v_old.raw();

        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                linear::Scalar uC = u(i, j);
                linear::Scalar uE = u(i + 1, j);
                linear::Scalar uW = u(i - 1, j);
                linear::Scalar uN = (j + 1 == Ny - 1) ? 1.0 : u(i, j + 1);  // lid u = 1
                linear::Scalar uS = u(i, j - 1);

                linear::Scalar conv_u = uC * (uE - uW) / (2.0 * dx) + v(i, j) * (uN - uS) / (2.0 * dy);
                linear::Scalar diff_u = (uE - 2.0 * uC + uW) / dx2 + (uN - 2.0 * uC + uS) / dy2;
                linear::Scalar dpdx = (p(i + 1, j) - p(i - 1, j)) / (2.0 * dx);

                ustar(i, j) = u_old(i, j) + dt * (-conv_u + nu * diff_u - dpdx);
            }
        }

        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                linear::Scalar vC = v(i, j);
                linear::Scalar vE = v(i + 1, j);
                linear::Scalar vW = v(i - 1, j);
                linear::Scalar vN = v(i, j + 1);
                linear::Scalar vS = v(i, j - 1);

                linear::Scalar conv_v = u(i, j) * (vE - vW) / (2.0 * dx) + vC * (vN - vS) / (2.0 * dy);
                linear::Scalar diff_v = (vE - 2.0 * vC + vW) / dx2 + (vN - 2.0 * vC + vS) / dy2;
                linear::Scalar dpdy = (p(i, j + 1) - p(i, j - 1)) / (2.0 * dy);

                vstar(i, j) = v_old(i, j) + dt * (-conv_v + nu * diff_v - dpdy);
            }
        }
    }

    // enforce lid-driven cavity BCs on velocity fields
    static void applyVelocityBC(core::Field &u, core::Field &v, int Nx, int Ny)
    {
        // u velocity: top lid = 1, others = 0
        for (int j = 0; j < Ny; ++j) {
            u(0,j) = 0.0;
            u(Nx-1,j) = 0.0;
        }
        for (int i = 0; i < Nx; ++i) {
            u(i,0) = 0.0;
            u(i,Ny-1) = 1.0;
        }
        // v velocity: all walls zero
        for (int j = 0; j < Ny; ++j) {
            v(0,j) = 0.0;
            v(Nx-1,j) = 0.0;
        }
        for (int i = 0; i < Nx; ++i) {
            v(i,0) = 0.0;
            v(i,Ny-1) = 0.0;
        }
    }

    static void solve(
        core::Mesh& mesh,
        core::Field& u,
        core::Field& v,
        core::Field& p,
        linear::Scalar nu,
        int maxIter,
        linear::Scalar alpha_u = 0.7,
        linear::Scalar alpha_v = 0.7,
        linear::Scalar alpha_p = 0.3,
        linear::Scalar time_step = 0.1)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        core::Field ustar = u;
        core::Field vstar = v;
        core::Field pCorr(mesh);

        applyVelocityBC(u, v, Nx, Ny);

        for (int iter = 0; iter < maxIter; ++iter) {
            auto sysU = assembly::MomentumAssembler::assembleU(mesh, u, v, p, nu, time_step);
            auto sysV = assembly::MomentumAssembler::assembleV(mesh, u, v, p, nu, time_step);

            linear::Vector solU =
                linear::BiCGSTABSolver().solve(sysU.A, sysU.b);
            linear::Vector solV =
                linear::BiCGSTABSolver().solve(sysV.A, sysV.b);

            // copy solution
            ustar.raw() = solU;
            vstar.raw() = solV;

            applyVelocityBC(ustar, vstar, Nx, Ny);

            auto pSys = assembly::PressureCorrectionAssembler::assemble(mesh, ustar, vstar);
            applyPressureBC(pCorr, mesh);
            linear::Vector solP = linear::BiCGSTABSolver().solve(pSys.A, pSys.b);
            pCorr.raw() = solP;

            // update with relaxation over whole field, BCs reapplied afterwards
            u.raw() += alpha_u * (ustar.raw() - u.raw());
            v.raw() += alpha_v * (vstar.raw() - v.raw());
            p.raw() += alpha_p * pCorr.raw();

            // velocity correction from pressure
            linear::Vector du = linear::zeros(mesh.totalCells());
            linear::Vector dv = linear::zeros(mesh.totalCells());
            for (int j = 1; j < Ny-1; ++j) {
                for (int i = 1; i < Nx-1; ++i) {
                    int P = core::id(i,j,Nx);
                    int E = core::id(i+1,j,Nx);
                    int Nn = core::id(i,j+1,Nx);
                    du(P) = - (pCorr.raw()(E) - pCorr.raw()(P)) / dx;
                    dv(P) = - (pCorr.raw()(Nn) - pCorr.raw()(P)) / dy;
                }
            }
            u.raw() += du;
            v.raw() += dv;

            applyVelocityBC(u, v, Nx, Ny);
        }
    }

private:
    static void applyPressureBC(core::Field &pCorr, const core::Mesh &mesh)
    {
        int Nx = mesh.Nx();
        int Ny = mesh.Ny();
        for (int j = 0; j < Ny; ++j) {
            pCorr(0,j) = 0.0;
            pCorr(Nx-1,j) = 0.0;
        }
        for (int i = 0; i < Nx; ++i) {
            pCorr(i,0) = 0.0;
            pCorr(i,Ny-1) = 0.0;
        }
    }
};

}