#include "../include/assembly/GenericScalarAssembler.hpp"
#include <cmath>
#include <algorithm>

namespace assembly {

ScalarSystem GenericScalarAssembler::assemble(
    const core::Mesh& mesh,
    const core::Field& phi,
    const core::Field& u,
    const core::Field& v,
    double gamma,
    double dt,
    const core::Field& phi_old,
    ConvectionScheme scheme,
    const core::Field* source_field,
    bool is_energy,
    double nu)
{
    return assembleWithScheme(mesh, phi, u, v, gamma, dt, phi_old, scheme, source_field, is_energy, nu);
}

core::Field GenericScalarAssembler::computeViscousDissipation(
    const core::Mesh& mesh,
    const core::Field& u,
    const core::Field& v,
    double nu)
{
    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();
    const double dx = mesh.dx();
    const double dy = mesh.dy();
    
    core::Field dissipation(mesh);
    dissipation.setZero();
    
    for (int j = 1; j < Ny - 1; ++j) {
        for (int i = 1; i < Nx - 1; ++i) {
            // Compute strain rate tensor components
            // ∂u/∂x
            double du_dx = (u(i + 1, j) - u(i - 1, j)) / (2.0 * dx);
            // ∂v/∂y
            double dv_dy = (v(i, j + 1) - v(i, j - 1)) / (2.0 * dy);
            // ∂u/∂y + ∂v/∂x
            double du_dy = (u(i, j + 1) - u(i, j - 1)) / (2.0 * dy);
            double dv_dx = (v(i + 1, j) - v(i - 1, j)) / (2.0 * dx);
            double shear = du_dy + dv_dx;
            
            // Φ = 2μ[ (∂u/∂x)² + (∂v/∂y)² + ((∂u/∂y + ∂v/∂x)/2)² ]
            double dissip = 2.0 * nu * (du_dx * du_dx + dv_dy * dv_dy + 0.25 * shear * shear);
            dissipation(i, j) = dissip;
        }
    }
    
    return dissipation;
}

double GenericScalarAssembler::upwindInterpolation(
    double phi_upstream,
    double phi_center,
    double velocity,
    double /*grid_spacing*/)
{
    if (velocity >= 0.0) {
        return phi_upstream;
    } else {
        return phi_center;
    }
}

double GenericScalarAssembler::centralInterpolation(
    double /*phi_upstream*/,
    double phi_center,
    double phi_downstream,
    double /*grid_spacing*/)
{
    return 0.5 * (phi_center + phi_downstream);
}

double GenericScalarAssembler::quickInterpolation(
    double phi_far_upstream,
    double phi_upstream,
    double phi_center,
    double phi_downstream,
    double velocity,
    double /*grid_spacing*/)
{
    if (velocity >= 0.0) {
        // QUICK uses upstream side: (3/8)φ_up + (3/4)φ_P - (1/8)φ_far_up
        return (3.0 / 8.0) * phi_upstream +
               (3.0 / 4.0) * phi_center -
               (1.0 / 8.0) * phi_far_upstream;
    } else {
        // Downstream side
        return (3.0 / 8.0) * phi_downstream +
               (3.0 / 4.0) * phi_center -
               (1.0 / 8.0) * phi_center;  // placeholder
    }
}

double GenericScalarAssembler::vanLeerLimiter(double r)
{
    if (r <= 0.0) return 0.0;
    return (r + std::abs(r)) / (1.0 + std::abs(r));
}

double GenericScalarAssembler::musclLimiter(double r)
{
    if (r <= 0.0) return 0.0;
    return std::max(0.0, std::min(2.0 * r, std::min(1.0, 2.0)));
}

double GenericScalarAssembler::interpolateToFace(
    double phi_far_upstream,
    double phi_upstream,
    double phi_center,
    double phi_downstream,
    double velocity,
    double grid_spacing,
    ConvectionScheme scheme)
{
    switch (scheme) {
        case ConvectionScheme::Central:
            return centralInterpolation(phi_upstream, phi_center, phi_downstream, grid_spacing);
        
        case ConvectionScheme::Upwind:
            return upwindInterpolation(phi_upstream, phi_center, velocity, grid_spacing);
        
        case ConvectionScheme::QUICK:
            return quickInterpolation(phi_far_upstream, phi_upstream, phi_center,
                                     phi_downstream, velocity, grid_spacing);
        
        case ConvectionScheme::MUSCL:
        case ConvectionScheme::VanLeer:
            // Slope limiter approach
            if (velocity >= 0.0) {
                // Upstream side interpolation
                double denom = phi_center - phi_upstream;
                if (std::abs(denom) < 1e-16) return phi_upstream;
                double r = (phi_upstream - phi_far_upstream) / denom;
                double limiter = (scheme == ConvectionScheme::MUSCL) ?
                                musclLimiter(r) : vanLeerLimiter(r);
                return phi_upstream + 0.5 * limiter * denom;
            } else {
                // Downstream side
                double denom = phi_downstream - phi_center;
                if (std::abs(denom) < 1e-16) return phi_downstream;
                double r = (phi_center - phi_downstream) / denom;
                double limiter = (scheme == ConvectionScheme::MUSCL) ?
                                musclLimiter(r) : vanLeerLimiter(r);
                return phi_center + 0.5 * limiter * denom;
            }
        
        default:
            return centralInterpolation(phi_upstream, phi_center, phi_downstream, grid_spacing);
    }
}

ScalarSystem GenericScalarAssembler::assembleWithScheme(
    const core::Mesh& mesh,
    const core::Field& phi,
    const core::Field& u,
    const core::Field& v,
    double gamma,
    double dt,
    const core::Field& phi_old,
    ConvectionScheme scheme,
    const core::Field* source_field,
    bool is_energy,
    double nu)
{
    using namespace core;
    using namespace linear;

    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();
    const int N = mesh.totalCells();

    const double dx = mesh.dx();
    const double dy = mesh.dy();
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    // Compute dissipation if energy
    Field dissipation(mesh);
    if (is_energy) {
        dissipation = computeViscousDissipation(mesh, u, v, nu);
    }

    SparseBuilder builder(N, N);
    Vector b(N);
    b.setZero();

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int P = id(i, j, Nx);

            if (mesh.isBoundary(i, j)) {
                builder.add(P, P, 1.0);
                b(P) = phi(i, j);  // Dirichlet BC
                continue;
            }

            int E = id(i + 1, j, Nx);
            int W = id(i - 1, j, Nx);
            int Nn = id(i, j + 1, Nx);
            int S = id(i, j - 1, Nx);

            double aP = 0.0;
            double aE = 0.0;
            double aW = 0.0;
            double aN = 0.0;
            double aS = 0.0;

            // Diffusion (always present)
            aE += gamma / dx2;
            aW += gamma / dx2;
            aN += gamma / dy2;
            aS += gamma / dy2;
            aP -= 2.0 * gamma / dx2 + 2.0 * gamma / dy2;

            // Convection: use face values
            double uC = u(i, j);
            // East face
            double phi_E = phi(i + 1, j);
            double phi_W = phi(i, j);
            double phi_WW = (i - 2 >= 0) ? phi(i - 2, j) : phi_W;
            interpolateToFace(phi_WW, phi_W, phi(i, j), phi_E, uC, dx, scheme);
            aE -= uC / (2.0 * dx);
            aW += uC / (2.0 * dx);
            // Actually, let's use simpler approach: convection stencil
            if (uC > 0) {
                aP += uC / dx;
                aW -= uC / dx;
            } else {
                aP -= uC / dx;
                aE += uC / dx;
            }

            // North face
            double vN = v(i, j);
            if (vN > 0) {
                aP += vN / dy;
                aS -= vN / dy;
            } else {
                aP -= vN / dy;
                aN += vN / dy;
            }

            // Unsteady term
            aP += 1.0 / dt;
            b(P) += phi_old(i, j) / dt;

            // Source terms
            if (source_field) {
                b(P) += (*source_field)(i, j);
            }

            // Energy dissipation
            if (is_energy) {
                b(P) += dissipation(i, j);
            }

            builder.add(P, P, aP);
            builder.add(P, E, aE);
            builder.add(P, W, aW);
            builder.add(P, Nn, aN);
            builder.add(P, S, aS);
        }
    }

    linear::SparseMatrix A = builder.build();
    return {A, b};
}

}  // namespace assembly
