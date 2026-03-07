#pragma once

#include "../core/Mesh.hpp"
#include "../core/Field.hpp"
#include "../core/Indexing.hpp"
#include "../linear/Vector.hpp"
#include "../linear/Matrix.hpp"

namespace assembly {

struct MomentumSystem {
    linear::SparseMatrix A;
    linear::Vector b;
};

class MomentumAssembler {
public:

    // assemble x–momentum (u) equation with TRANSIENT unsteady term
    // properly using OLD velocity (u_old, v_old) from previous time step
    static MomentumSystem assembleUTransient(
        const core::Mesh& mesh,
        const core::Field& u,        // current iterate
        const core::Field& u_old,    // previous time step
        const core::Field& v,        // current iterate
        const core::Field& v_old,    // previous time step
        const core::Field& p,
        double nu,
        double dt)                   // true time step
    {
        using namespace core;
        using namespace linear;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int N  = mesh.totalCells();

        const double dx = mesh.dx();
        const double dy = mesh.dy();

        const double dx2 = dx * dx;
        const double dy2 = dy * dy;

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

                double uP = u(i,j);
                double vP = v(i,j);

                double aP = 0.0;
                double aE = 0.0;
                double aW = 0.0;
                double aN = 0.0;
                double aS = 0.0;

                // Diffusion
                aE += nu/dx2;
                aW += nu/dx2;
                aN += nu/dy2;
                aS += nu/dy2;
                aP -= 2.0*nu/dx2 + 2.0*nu/dy2;

                // Convection (central)
                aE -= uP/(2.0*dx);
                aW += uP/(2.0*dx);
                aN -= vP/(2.0*dy);
                aS += vP/(2.0*dy);

                // Unsteady term: (u - u_old) / dt
                aP += 1.0/dt;
                double u_old_val = u_old(i,j);
                b(P) += u_old_val / dt;

                // add a tiny constant source if field is exactly zero
                if (u_old_val == 0.0 && uP == 0.0) {
                    b(P) += 1e-6;  // small kick to break symmetry
                }

                builder.add(P, P, aP);
                builder.add(P, E, aE);
                builder.add(P, W, aW);
                // North: if lid (u=1), put aN*u_lid on RHS and do not add (P,Nn) to matrix
                if (mesh.isTopBoundary(j + 1)) {
                    b(P) += aN * 1.0;
                } else {
                    builder.add(P, Nn, aN);
                }
                builder.add(P, S, aS);

                // pressure gradient term - to RHS
                b(P) -= (p(i+1,j) - p(i-1,j)) / (2.0*dx);
            }
        }

        SparseMatrix A = builder.build();
        return {A, b};
    }

    // assemble y–momentum (v) equation with TRANSIENT unsteady term
    static MomentumSystem assembleVTransient(
        const core::Mesh& mesh,
        const core::Field& u,        // current iterate
        const core::Field& u_old,    // previous time step
        const core::Field& v,        // current iterate
        const core::Field& v_old,    // previous time step
        const core::Field& p,
        double nu,
        double dt)                   // true time step
    {
        using namespace core;
        using namespace linear;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int N  = mesh.totalCells();

        const double dx = mesh.dx();
        const double dy = mesh.dy();

        const double dx2 = dx * dx;
        const double dy2 = dy * dy;

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
                    b(P) = 0.0;  // v=0 at all boundaries
                    continue;
                }

                int E = id(i+1, j, Nx);
                int W = id(i-1, j, Nx);
                int Nn = id(i, j+1, Nx);
                int S = id(i, j-1, Nx);

                double uP = u(i,j);
                double vP = v(i,j);

                double aP = 0.0;
                double aE = 0.0;
                double aW = 0.0;
                double aN = 0.0;
                double aS = 0.0;

                // Diffusion
                aE += nu/dx2;
                aW += nu/dx2;
                aN += nu/dy2;
                aS += nu/dy2;
                aP -= 2.0*nu/dx2 + 2.0*nu/dy2;

                // Convection (central)
                aE -= uP/(2.0*dx);
                aW += uP/(2.0*dx);
                aN -= vP/(2.0*dy);
                aS += vP/(2.0*dy);

                // Unsteady term: (v - v_old) / dt
                aP += 1.0/dt;
                double v_old_val = v_old(i,j);
                b(P) += v_old_val / dt;

                // small kick if both old and current are zero
                if (v_old_val == 0.0 && vP == 0.0) {
                    b(P) += 1e-6;
                }

                builder.add(P, P, aP);
                builder.add(P, E, aE);
                builder.add(P, W, aW);
                builder.add(P, Nn, aN);
                builder.add(P, S, aS);

                // pressure gradient term (y direction)
                b(P) -= (p(i,j+1) - p(i,j-1)) / (2.0*dy);
            }
        }

        SparseMatrix A = builder.build();
        return {A, b};
    }

    // assemble x–momentum (u) equation with optional pressure gradient
    static MomentumSystem assembleU(
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v,
        const core::Field& p,
        double nu,
        double time_step = 0.0)  // pseudo-time step for stabilization
    {
        using namespace core;
        using namespace linear;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int N  = mesh.totalCells();

        const double dx = mesh.dx();
        const double dy = mesh.dy();

        const double dx2 = dx * dx;
        const double dy2 = dy * dy;

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

                double uP = u(i,j);
                double vP = v(i,j);

                double aP = 0.0;
                double aE = 0.0;
                double aW = 0.0;
                double aN = 0.0;
                double aS = 0.0;

                // Diffusion
                aE += nu/dx2;
                aW += nu/dx2;
                aN += nu/dy2;
                aS += nu/dy2;
                aP -= 2.0*nu/dx2 + 2.0*nu/dy2;

                // Convection (central)
                aE -= uP/(2.0*dx);
                aW += uP/(2.0*dx);
                aN -= vP/(2.0*dy);
                aS += vP/(2.0*dy);

                builder.add(P, P, aP);
                builder.add(P, E, aE);
                builder.add(P, W, aW);
                builder.add(P, Nn, aN);
                builder.add(P, S, aS);

                // pressure gradient term - to RHS
                b(P) = - (p(i+1,j) - p(i-1,j)) / (2.0*dx);

                // Unsteady term (pseudo-time for stabilization)
                if (time_step > 0.0) {
                    aP += 1.0/time_step;
                    b(P) += uP/time_step;
                }
            }
        }

        SparseMatrix A = builder.build();
        return {A, b};
    }

    // assemble y–momentum (v) equation with pressure gradient
    static MomentumSystem assembleV(
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v,
        const core::Field& p,
        double nu,
        double time_step = 0.0)  // pseudo-time step for stabilization
    {
        using namespace core;
        using namespace linear;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int N  = mesh.totalCells();

        const double dx = mesh.dx();
        const double dy = mesh.dy();

        const double dx2 = dx * dx;
        const double dy2 = dy * dy;

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
                    b(P) = 0.0;  // v=0 at all boundaries
                    continue;
                }

                int E = id(i+1, j, Nx);
                int W = id(i-1, j, Nx);
                int Nn = id(i, j+1, Nx);
                int S = id(i, j-1, Nx);

                double uP = u(i,j);
                double vP = v(i,j);

                double aP = 0.0;
                double aE = 0.0;
                double aW = 0.0;
                double aN = 0.0;
                double aS = 0.0;

                // Diffusion
                aE += nu/dx2;
                aW += nu/dx2;
                aN += nu/dy2;
                aS += nu/dy2;
                aP -= 2.0*nu/dx2 + 2.0*nu/dy2;

                // Convection (central)
                aE -= uP/(2.0*dx);
                aW += uP/(2.0*dx);
                aN -= vP/(2.0*dy);
                aS += vP/(2.0*dy);

                builder.add(P, P, aP);
                builder.add(P, E, aE);
                builder.add(P, W, aW);
                builder.add(P, Nn, aN);
                builder.add(P, S, aS);

                // pressure gradient term (y direction)
                b(P) = - (p(i,j+1) - p(i,j-1)) / (2.0*dy);

                // Unsteady term (pseudo-time for stabilization)
                if (time_step > 0.0) {
                    aP += 1.0/time_step;
                    b(P) += vP/time_step;
                }
            }
        }

        SparseMatrix A = builder.build();
        return {A, b};
    }
};

}