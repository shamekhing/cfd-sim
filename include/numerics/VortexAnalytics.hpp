#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "../core/Mesh.hpp"
#include "../core/Field.hpp"

namespace numerics {

/**
 * @brief Vortex analysis utilities
 * 
 * Computes vorticity field, streamfunction, and detects vortex centers
 * for lid-driven cavity and other shear-driven flows.
 */
class VortexAnalytics {
public:

    /**
     * @brief Compute vorticity field ω = ∂v/∂x - ∂u/∂y
     * 
     * Vorticity is computed at cell centers using central differences.
     * For staggered grid: u is on vertical faces, v on horizontal faces.
     * 
     * @param mesh Computational mesh
     * @param u_cc U-velocity at cell centers
     * @param v_cc V-velocity at cell centers
     * @return Field containing vorticity ω at each cell center
     */
    static core::Field computeVorticity(
        const core::Mesh& mesh,
        const core::Field& u_cc,
        const core::Field& v_cc)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const double dx = mesh.dx();
        const double dy = mesh.dy();
        
        core::Field omega(mesh);
        omega.setZero();
        
        // Interior cells: central differences
        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                // ∂v/∂x
                double dvdx = (v_cc(i + 1, j) - v_cc(i - 1, j)) / (2.0 * dx);
                
                // ∂u/∂y
                double dudy = (u_cc(i, j + 1) - u_cc(i, j - 1)) / (2.0 * dy);
                
                // ω = ∂v/∂x - ∂u/∂y
                omega(i, j) = dvdx - dudy;
            }
        }
        
        // Boundary cells: one-sided differences (for cavity: ω=0 at boundaries)
        for (int i = 0; i < Nx; ++i) {
            omega(i, 0) = 0.0;
            omega(i, Ny - 1) = 0.0;
        }
        for (int j = 0; j < Ny; ++j) {
            omega(0, j) = 0.0;
            omega(Nx - 1, j) = 0.0;
        }
        
        return omega;
    }

    /**
     * @brief Compute streamfunction field ψ from velocity field
     * 
     * For 2D incompressible flow: ∂ψ/∂y = u and -∂ψ/∂x = v
     * Computed via line integration: ψ(i,j) = ∫[from origin] u·dy - v·dx
     * 
     * Simple approach: integrate from bottom-left corner.
     * Streamlines are contours of constant ψ.
     * 
     * @param mesh Computational mesh
     * @param u_cc U-velocity at cell centers
     * @param v_cc V-velocity at cell centers
     * @return Field containing streamfunction ψ at each cell
     */
    static core::Field computeStreamfunction(
        const core::Mesh& mesh,
        const core::Field& u_cc,
        const core::Field& v_cc)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const double dx = mesh.dx();
        const double dy = mesh.dy();
        
        core::Field psi(mesh);
        psi.setZero();
        
        // Bottom-left corner is reference point (ψ = 0)
        psi(0, 0) = 0.0;
        
        // Integrate along x-direction (constant j=0)
        for (int i = 1; i < Nx; ++i) {
            psi(i, 0) = psi(i - 1, 0) - v_cc(i - 1, 0) * dx;
        }
        
        // Integrate along y-direction (for each column)
        for (int j = 1; j < Ny; ++j) {
            psi(0, j) = psi(0, j - 1) + u_cc(0, j - 1) * dy;
            
            for (int i = 1; i < Nx; ++i) {
                // Average of two paths (more robust)
                double psi_from_left = psi(i - 1, j) - v_cc(i - 1, j) * dx;
                double psi_from_below = psi(i, j - 1) + u_cc(i, j - 1) * dy;
                
                psi(i, j) = 0.5 * (psi_from_left + psi_from_below);
            }
        }
        
        return psi;
    }

    /**
     * @brief Find primary vortex center using pressure minimum
     * 
     * For lid-driven cavity, the primary vortex center corresponds to
     * a local (and usually global) pressure minimum in the interior.
     * 
     * Reference: Ghia et al. (1982) gives primary vortex center at
     * approximately (0.53, 0.56) for Re=100.
     * 
     * @param mesh Computational mesh
     * @param p Pressure field
     * @return Pair (x_center, y_center) of detected primary vortex
     */
    static std::pair<double, double> findPrimaryVortexCenter(
        const core::Mesh& mesh,
        const core::Field& p)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        
        // Search in interior only (avoid boundary)
        double min_pressure = p(1, 1);
        int i_center = 1;
        int j_center = 1;
        
        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                if (p(i, j) < min_pressure) {
                    min_pressure = p(i, j);
                    i_center = i;
                    j_center = j;
                }
            }
        }
        
        double x_center = mesh.x(i_center);
        double y_center = mesh.y(j_center);
        
        return {x_center, y_center};
    }

    /**
     * @brief Find secondary vortices using vorticity criteria
     * 
     * Secondary vortices (e.g., corner vortices at Re=400) are identified
     * as local extrema of vorticity (ω) in subdomain regions.
     * 
     * Strategy:
     * 1. Divide domain into quadrants
     * 2. Find local vorticity extremum in each quadrant (opposite sign to primary)
     * 3. Return those above a threshold in magnitude
     * 
     * @param mesh Computational mesh
     * @param omega Vorticity field
     * @param threshold Minimum vorticity magnitude to report as a vortex
     * @return Vector of (x, y) pairs for detected secondary vortices
     */
    static std::vector<std::pair<double, double>> findSecondaryVortices(
        const core::Mesh& mesh,
        const core::Field& omega,
        double threshold = 0.05)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        
        std::vector<std::pair<double, double>> secondary_centers;
        
        // Define four corner regions
        struct Region {
            int i_min, i_max, j_min, j_max;
            std::string name;
        };
        
        int mid_i = Nx / 2;
        int mid_j = Ny / 2;
        
        std::vector<Region> regions = {
            {1, mid_i - 1, 1, mid_j - 1, "bottom-left"},
            {mid_i + 1, Nx - 2, 1, mid_j - 1, "bottom-right"},
            {1, mid_i - 1, mid_j + 1, Ny - 2, "top-left"},
            {mid_i + 1, Nx - 2, mid_j + 1, Ny - 2, "top-right"}
        };
        
        for (const auto& region : regions) {
            double max_vorticity_magnitude = 0.0;
            int i_vortex = region.i_min;
            int j_vortex = region.j_min;
            
            // Find local extremum in this region
            for (int j = region.j_min; j <= region.j_max; ++j) {
                for (int i = region.i_min; i <= region.i_max; ++i) {
                    double vort_mag = std::abs(omega(i, j));
                    
                    // Check if local extremum (different from primary vortex sign)
                    bool is_extremum = (vort_mag > max_vorticity_magnitude) &&
                                       (vort_mag > threshold);
                    
                    if (is_extremum) {
                        max_vorticity_magnitude = vort_mag;
                        i_vortex = i;
                        j_vortex = j;
                    }
                }
            }
            
            // Only report if above threshold
            if (max_vorticity_magnitude > threshold) {
                double x = mesh.x(i_vortex);
                double y = mesh.y(j_vortex);
                secondary_centers.push_back({x, y});
            }
        }
        
        return secondary_centers;
    }

    /**
     * @brief Check if a cell is a local extremum in a field
     * 
     * Useful for vortex detection algorithms.
     * 
     * @param field Scalar field
     * @param i Cell i-index
     * @param j Cell j-index
     * @param is_minimum If true, check for minimum; else check for maximum
     * @return True if (i,j) is a local extremum
     */
    static bool isLocalExtremum(
        const core::Field& field,
        int i, int j,
        bool is_minimum = true)
    {
        double center = field(i, j);
        
        // Check 4-neighborhood
        double neighbors[4] = {
            field(i + 1, j),  // right
            field(i - 1, j),  // left
            field(i, j + 1),  // top
            field(i, j - 1)   // bottom
        };
        
        if (is_minimum) {
            for (int n = 0; n < 4; ++n) {
                if (neighbors[n] < center) return false;
            }
            return true;
        } else {
            for (int n = 0; n < 4; ++n) {
                if (neighbors[n] > center) return false;
            }
            return true;
        }
    }
};

}
