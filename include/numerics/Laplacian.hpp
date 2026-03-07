#pragma once

#include "../linear/Vector.hpp"
#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"

namespace numerics {

/**
 * @brief 5-point Laplacian operator on cell-centered grid
 * 
 * Computes ∇²φ = ∂²φ/∂x² + ∂²φ/∂y² using central differences:
 * 
 * ∇²φ(i,j) ≈ [φ(i+1,j) - 2φ(i,j) + φ(i-1,j)] / (Δx)² 
 *           + [φ(i,j+1) - 2φ(i,j) + φ(i,j-1)] / (Δy)²
 * 
 * At boundaries: uses one-sided differences or Dirichlet values
 */
class Laplacian {
public:
    /**
     * @brief Compute Laplacian at a single cell (i,j)
     * @return ∇²φ at cell center
     */
    static linear::Scalar operator5Point(
        const core::Field& phi,
        const core::Mesh& mesh,
        int i, int j)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        
        // Clamp to interior to avoid out-of-bounds
        i = std::max(0, std::min(i, Nx - 1));
        j = std::max(0, std::min(j, Ny - 1));
        
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();
        const linear::Scalar dx2 = dx * dx;
        const linear::Scalar dy2 = dy * dy;
        
        linear::Scalar laplace = 0.0;
        
        // X-direction: (φ_E - 2φ_P + φ_W) / dx²
        if (i > 0 && i < Nx - 1) {
            laplace += (phi(i + 1, j) - 2.0 * phi(i, j) + phi(i - 1, j)) / dx2;
        } else if (i == 0) {
            // One-sided: (φ_center - φ_W) as proxy
            laplace += (phi(i + 1, j) - phi(i, j)) / dx2;
        } else if (i == Nx - 1) {
            // One-sided: (φ_E - φ_center) as proxy
            laplace += (phi(i, j) - phi(i - 1, j)) / dx2;
        }
        
        // Y-direction: (φ_N - 2φ_P + φ_S) / dy²
        if (j > 0 && j < Ny - 1) {
            laplace += (phi(i, j + 1) - 2.0 * phi(i, j) + phi(i, j - 1)) / dy2;
        } else if (j == 0) {
            // One-sided
            laplace += (phi(i, j + 1) - phi(i, j)) / dy2;
        } else if (j == Ny - 1) {
            // One-sided
            laplace += (phi(i, j) - phi(i, j - 1)) / dy2;
        }
        
        return laplace;
    }
    
    /**
     * @brief Compute Laplacian field: ∇²φ everywhere
     * @return Field containing Laplacian values at each cell
     */
    static core::Field apply(
        const core::Field& phi,
        const core::Mesh& mesh)
    {
        core::Field laplace(mesh);
        
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                laplace(i, j) = operator5Point(phi, mesh, i, j);
            }
        }
        
        return laplace;
    }
};

}  // namespace numerics
