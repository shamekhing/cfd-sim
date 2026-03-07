#pragma once

#include "../linear/Vector.hpp"
#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"
#include "../../include/core/UField.hpp"
#include "../../include/core/VField.hpp"

namespace numerics {

/**
 * @brief Strain rate tensor utilities for turbulence modeling
 * 
 * For 2D incompressible flow, the strain rate tensor is:
 *   S_ij = 0.5 * (∂u_i/∂x_j + ∂u_j/∂x_i)
 * 
 * Components:
 *   S_xx = ∂u/∂x
 *   S_yy = ∂v/∂y
 *   S_xy = S_yx = 0.5 * (∂u/∂y + ∂v/∂x)
 * 
 * Magnitude squared (Frobenius norm):
 *   ||S||² = S_xx² + S_yy² + 2*S_xy²
 * 
 * For turbulence production:
 *   P_k = 2*μ_t*||S||²
 */
class StrainRate
{
public:
    /**
     * @brief Compute strain rate magnitude squared at cell centers
     * 
     * Returns: Field containing 2*||S||² = 2*(S_xx² + S_yy² + 2*S_xy²)
     * This is the term used in production: P_k = μ_t * (2*||S||²)
     * 
     * @param u U-velocity field (staggered, on vertical faces)
     * @param v V-velocity field (staggered, on horizontal faces)
     * @param mesh Computational mesh
     * @return Field with 2*||S||² at each cell center
     */
    static core::Field computeMagnitudeSquared(
        const core::UField& u,
        const core::VField& v,
        const core::Mesh& mesh)
    {
        using namespace core;
        
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();
        
        Field result(mesh);
        
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                // Compute velocity gradients at cell center (i, j)
                
                // ∂u/∂x at cell center: central difference using face velocities
                // u is stored at vertical faces: u(i,j) is west face, u(i+1,j) is east face
                linear::Scalar dudx = (u(i + 1, j) - u(i, j)) / dx;
                
                // ∂v/∂y at cell center: central difference using face velocities
                // v is stored at horizontal faces: v(i,j) is south face, v(i,j+1) is north face
                linear::Scalar dvdy = (v(i, j + 1) - v(i, j)) / dy;
                
                // ∂u/∂y at cell center: requires averaging u to cell centers first, then difference
                // Average u to north and south cell centers, then central difference
                linear::Scalar u_center_south = linear::Scalar(0.5) * (u(i, j) + u(i + 1, j));
                linear::Scalar u_center_north = (j + 1 < Ny) ? linear::Scalar(0.5) * (u(i, j + 1) + u(i + 1, j + 1)) : u_center_south;
                linear::Scalar dudy = (u_center_north - u_center_south) / dy;
                
                // ∂v/∂x at cell center: requires averaging v to cell centers first, then difference
                // Average v to west and east cell centers, then central difference
                linear::Scalar v_center_west = linear::Scalar(0.5) * (v(i, j) + v(i, j + 1));
                linear::Scalar v_center_east = (i + 1 < Nx) ? linear::Scalar(0.5) * (v(i + 1, j) + v(i + 1, j + 1)) : v_center_west;
                linear::Scalar dvdx = (v_center_east - v_center_west) / dx;
                
                // Strain rate components
                linear::Scalar S_xx = dudx;
                linear::Scalar S_yy = dvdy;
                linear::Scalar S_xy = linear::Scalar(0.5) * (dudy + dvdx);
                
                // Magnitude squared: ||S||² = S_xx² + S_yy² + 2*S_xy²
                // Return 2*||S||² for production term: P_k = μ_t * (2*||S||²)
                linear::Scalar magnitude_sq = linear::Scalar(2.0) * (S_xx * S_xx + S_yy * S_yy + linear::Scalar(2.0) * S_xy * S_xy);
                
                result(i, j) = magnitude_sq;
            }
        }
        
        return result;
    }
    
    /**
     * @brief Compute individual strain rate components at cell centers
     * 
     * @param u U-velocity field
     * @param v V-velocity field
     * @param mesh Computational mesh
     * @param[out] S_xx ∂u/∂x component
     * @param[out] S_yy ∂v/∂y component
     * @param[out] S_xy 0.5*(∂u/∂y + ∂v/∂x) component
     */
    static void computeComponents(
        const core::UField& u,
        const core::VField& v,
        const core::Mesh& mesh,
        core::Field& S_xx,
        core::Field& S_yy,
        core::Field& S_xy)
    {
        using namespace core;
        
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();
        
        S_xx.resize(mesh);
        S_yy.resize(mesh);
        S_xy.resize(mesh);
        
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                // ∂u/∂x
                linear::Scalar dudx = (u(i + 1, j) - u(i, j)) / dx;
                
                // ∂v/∂y
                linear::Scalar dvdy = (v(i, j + 1) - v(i, j)) / dy;
                
                // ∂u/∂y
                linear::Scalar u_center_south = linear::Scalar(0.5) * (u(i, j) + u(i + 1, j));
                linear::Scalar u_center_north = (j + 1 < Ny) ? linear::Scalar(0.5) * (u(i, j + 1) + u(i + 1, j + 1)) : u_center_south;
                linear::Scalar dudy = (u_center_north - u_center_south) / dy;
                
                // ∂v/∂x
                linear::Scalar v_center_west = linear::Scalar(0.5) * (v(i, j) + v(i, j + 1));
                linear::Scalar v_center_east = (i + 1 < Nx) ? linear::Scalar(0.5) * (v(i + 1, j) + v(i + 1, j + 1)) : v_center_west;
                linear::Scalar dvdx = (v_center_east - v_center_west) / dx;
                
                S_xx(i, j) = dudx;
                S_yy(i, j) = dvdy;
                S_xy(i, j) = linear::Scalar(0.5) * (dudy + dvdx);
            }
        }
    }
};

}  // namespace numerics
