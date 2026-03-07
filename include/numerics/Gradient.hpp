#pragma once

#include "../linear/Vector.hpp"
#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"
#include "../../include/core/UField.hpp"
#include "../../include/core/VField.hpp"

namespace numerics {

/**
 * @brief Gradient operator for scalar fields
 * 
 * Computes ∇φ = (∂φ/∂x, ∂φ/∂y) using central differences
 */
class Gradient {
public:
    /**
     * @brief Compute ∂φ/∂x at a location (i, j)
     * Uses central differences: ∂φ/∂x ≈ (φ(i+1,j) - φ(i-1,j)) / (2Δx)
     */
    static linear::Scalar dphidx(
        const core::Field& phi,
        const core::Mesh& mesh,
        int i, int j)
    {
        const int Nx = mesh.Nx();
        const linear::Scalar dx = mesh.dx();
        
        // Clamp i to [0, Nx-1]
        i = std::max(0, std::min(i, Nx - 1));
        
        if (i == 0) {
            // One-sided forward: (φ_1 - φ_0) / dx
            return (phi(1, j) - phi(0, j)) / dx;
        } else if (i == Nx - 1) {
            // One-sided backward: (φ_N-1 - φ_N-2) / dx
            return (phi(Nx - 1, j) - phi(Nx - 2, j)) / dx;
        } else {
            // Central difference
            return (phi(i + 1, j) - phi(i - 1, j)) / (2.0 * dx);
        }
    }
    
    /**
     * @brief Compute ∂φ/∂y at a location (i, j)
     * Uses central differences: ∂φ/∂y ≈ (φ(i,j+1) - φ(i,j-1)) / (2Δy)
     */
    static linear::Scalar dphidy(
        const core::Field& phi,
        const core::Mesh& mesh,
        int i, int j)
    {
        const int Ny = mesh.Ny();
        const linear::Scalar dy = mesh.dy();
        
        // Clamp j to [0, Ny-1]
        j = std::max(0, std::min(j, Ny - 1));
        
        if (j == 0) {
            // One-sided forward: (φ_1 - φ_0) / dy
            return (phi(i, 1) - phi(i, 0)) / dy;
        } else if (j == Ny - 1) {
            // One-sided backward: (φ_N-1 - φ_N-2) / dy
            return (phi(i, Ny - 1) - phi(i, Ny - 2)) / dy;
        } else {
            // Central difference
            return (phi(i, j + 1) - phi(i, j - 1)) / (2.0 * dy);
        }
    }
    
    /**
     * @brief Compute ∇φ field: ∂φ/∂x array at each cell center
     */
    static core::Field gradientX(
        const core::Field& phi,
        const core::Mesh& mesh)
    {
        core::Field dphi_dx(mesh);
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                dphi_dx(i, j) = dphidx(phi, mesh, i, j);
            }
        }
        
        return dphi_dx;
    }
    
    /**
     * @brief Compute ∇φ field: ∂φ/∂y array at each cell center
     */
    static core::Field gradientY(
        const core::Field& phi,
        const core::Mesh& mesh)
    {
        core::Field dphi_dy(mesh);
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                dphi_dy(i, j) = dphidy(phi, mesh, i, j);
            }
        }
        
        return dphi_dy;
    }
    
    /**
     * @brief Compute pressure gradient at a U-face (i, j) in x-direction
     * Pressure gradient ∂p/∂x at east face of cell (i-1, j)
     * Interpolates: (p_E - p_W) / (2Δx)
     */
    static linear::Scalar pressureGradientAtUFace(
        const core::Field& p,
        const core::Mesh& mesh,
        int i, int j)
    {
        const linear::Scalar dx = mesh.dx();
        
        // Face (i,j) is the east face of cell (i-1, j)
        // Pressure is at cell centers, so use cells (i-1,j) and (i,j)
        int iL = i - 1;  // Cell to the left
        int iR = i;      // Cell to the right
        
        return (p(iR, j) - p(iL, j)) / dx;
    }
    
    /**
     * @brief Compute pressure gradient at a V-face (i, j) in y-direction
     * Pressure gradient ∂p/∂y at north face of cell (i, j-1)
     * Interpolates: (p_N - p_S) / (2Δy)
     */
    static linear::Scalar pressureGradientAtVFace(
        const core::Field& p,
        const core::Mesh& mesh,
        int i, int j)
    {
        const linear::Scalar dy = mesh.dy();
        
        // Face (i,j) is the north face of cell (i, j-1)
        // Pressure is at cell centers, so use cells (i,j-1) and (i,j)
        int jB = j - 1;  // Cell below
        int jT = j;      // Cell above
        
        return (p(i, jT) - p(i, jB)) / dy;
    }
};

}  // namespace numerics
