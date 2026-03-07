#pragma once

#include "../linear/Vector.hpp"
#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"
#include "../../include/core/UField.hpp"
#include "../../include/core/VField.hpp"

namespace numerics {

/**
 * @brief Divergence operator for velocity fields
 * 
 * Computes ∇·u = ∂u/∂x + ∂v/∂y at cell centers
 * Used for:
 * - Mass conservation check: ∇·u ≈ 0 for incompressible flow
 * - Pressure correction assembly
 */
class Divergence {
public:
    /**
     * @brief Compute divergence ∇·u at cell center (i, j)
     * 
     * On staggered grid:
     * - U-velocity u is at vertical faces (east/west of cell)
     * - V-velocity v is at horizontal faces (north/south of cell)
     * 
     * Flux balance over cell (i,j):
     * ∇·u ≈ [u(i+1,j) - u(i,j)] / Δx + [v(i,j+1) - v(i,j)] / Δy
     * 
     * This is a volume-integral form: (outflow - inflow) / volume
     */
    static linear::Scalar divergenceAtCell(
        const core::UField& u,
        const core::VField& v,
        const core::Mesh& mesh,
        int i, int j)
    {
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();
        
        // U-flux: u enters from west face (i, j), exits from east face (i+1, j)
        linear::Scalar div_u = (u(i + 1, j) - u(i, j)) / dx;
        
        // V-flux: v enters from south face (i, j), exits from north face (i, j+1)
        linear::Scalar div_v = (v(i, j + 1) - v(i, j)) / dy;
        
        return div_u + div_v;
    }
    
    /**
     * @brief Compute divergence field: ∇·u at all cell centers
     * @return Field with divergence values (should be ≈ 0 for incompressible)
     */
    static core::Field apply(
        const core::UField& u,
        const core::VField& v,
        const core::Mesh& mesh)
    {
        core::Field div(mesh);
        
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                div(i, j) = divergenceAtCell(u, v, mesh, i, j);
            }
        }
        
        return div;
    }
    
    /**
     * @brief Compute L2 norm of divergence field
     * Measure of how well mass is conserved: should approach 0
     */
    static linear::Scalar divergenceL2Norm(
        const core::UField& u,
        const core::VField& v,
        const core::Mesh& mesh)
    {
        core::Field div = apply(u, v, mesh);
        return div.norm();
    }
    
    /**
     * @brief Compute max absolute divergence
     * Alternative convergence measure
     */
    static linear::Scalar divergenceMaxAbs(
        const core::UField& u,
        const core::VField& v,
        const core::Mesh& mesh)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        
        linear::Scalar max_div = 0.0;
        
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                linear::Scalar div = divergenceAtCell(u, v, mesh, i, j);
                max_div = std::max(max_div, std::abs(div));
            }
        }
        
        return max_div;
    }
};

}  // namespace numerics
