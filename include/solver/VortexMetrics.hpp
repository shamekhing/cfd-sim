#pragma once

#include <cmath>
#include <vector>
#include <utility>
#include "../core/Mesh.hpp"
#include "../core/Field.hpp"
#include "../numerics/VortexAnalytics.hpp"

namespace solver {

/**
 * @brief Vortex validation result structure
 */
struct VortexValidationResult {
    double center_x = 0.0;
    double center_y = 0.0;
    double circulation = 0.0;
    double enstrophy = 0.0;
    double max_vorticity = 0.0;
    
    // Error metrics (vs reference)
    double center_error = 0.0;        // L2 distance to reference center
    double circulation_error = 0.0;   // Relative % error
};

/**
 * @brief Wall shear profile result
 */
struct WallShearProfile {
    std::vector<double> x_coords;      // Wall-aligned coordinate
    std::vector<double> tau_w;         // Shear stress magnitude
    double tau_w_max = 0.0;
    double tau_w_mean = 0.0;
};

/**
 * @brief Vortex metrics computation and validation
 * 
 * Extends ResidualMonitor pattern to compute vortex-related quantities
 * and validate against reference data.
 */
class VortexMetrics {
public:

    /**
     * @brief Compute circulation around a vortex core
     * 
     * Γ = ∮ ω dA over circular contour around vortex center
     * 
     * For 2D: integrated vorticity (ω) over a disk of given radius
     * 
     * @param mesh Computational mesh
     * @param omega Vorticity field
     * @param center_x X-coordinate of vortex center
     * @param center_y Y-coordinate of vortex center
     * @param radius Radius of integration contour
     * @return Circulation magnitude
     */
    static double computeCirculation(
        const core::Mesh& mesh,
        const core::Field& omega,
        double center_x,
        double center_y,
        double radius = 0.1)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const double dx = mesh.dx();
        const double dy = mesh.dy();
        
        double circulation = 0.0;
        int count = 0;
        
        // Integrate vorticity over disk
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = mesh.x(i);
                double y = mesh.y(j);
                
                double dist = std::sqrt((x - center_x) * (x - center_x) +
                                       (y - center_y) * (y - center_y));
                
                if (dist <= radius) {
                    circulation += std::abs(omega(i, j)) * dx * dy;
                    count++;
                }
            }
        }
        
        return circulation;
    }

    /**
     * @brief Compute enstrophy (integrated vorticity squared)
     * 
     * E = 0.5 ∫ ω² dA
     * 
     * Global measure of rotational kinetic energy.
     * 
     * @param mesh Computational mesh
     * @param omega Vorticity field
     * @return Enstrophy magnitude
     */
    static double computeEnstrophy(
        const core::Mesh& mesh,
        const core::Field& omega)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const double dx = mesh.dx();
        const double dy = mesh.dy();
        
        double enstrophy = 0.0;
        
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double omega_sq = omega(i, j) * omega(i, j);
                enstrophy += 0.5 * omega_sq * dx * dy;
            }
        }
        
        return enstrophy;
    }

    /**
     * @brief Compute max vorticity magnitude in field
     * 
     * @param omega Vorticity field
     * @return Maximum |ω| value
     */
    static double computeMaxVorticity(const core::Field& omega) {
        double max_omega = 0.0;
        
        for (int idx = 0; idx < omega.size(); ++idx) {
            max_omega = std::max(max_omega, std::abs(omega[idx]));
        }
        
        return max_omega;
    }

    /**
     * @brief Compute wall shear stress profile along a wall
     * 
     * τ_w = μ (∂u/∂n) where n is normal to wall
     * 
     * For lid-driven cavity:
     * - Top wall (y=1): τ_w = ν ∂u/∂y
     * - Bottom wall (y=0): τ_w = -ν ∂u/∂y
     * - Left/right: τ_w = ν ∂v/∂x
     * 
     * @param mesh Computational mesh
     * @param u_cc U-velocity
     * @param v_cc V-velocity
     * @param nu Kinematic viscosity
     * @param nu_t Turbulent viscosity (optional, default 0)
     * @return WallShearProfile
     */
    static WallShearProfile computeWallShearProfile(
        const core::Mesh& mesh,
        const core::Field& u_cc,
        const core::Field& v_cc [[maybe_unused]],
        double nu,
        const core::Field* nu_t = nullptr)
    {
        WallShearProfile profile;
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const double dy = mesh.dy();
        
        // Top wall (y = 1): u-velocity gradient
        profile.x_coords.clear();
        profile.tau_w.clear();
        
        int j_top = Ny - 1;
        int j_near_top = Ny - 2;
        
        for (int i = 0; i < Nx; ++i) {
            double x = mesh.x(i);
            
            // One-sided difference at top wall
            double dudy = (u_cc(i, j_top) - u_cc(i, j_near_top)) / dy;
            
            double nu_eff = nu;
            if (nu_t != nullptr) {
                nu_eff += (*nu_t)(i, j_top);
            }
            
            double tau = nu_eff * dudy;
            
            profile.x_coords.push_back(x);
            profile.tau_w.push_back(std::abs(tau));
        }
        
        // Statistics
        profile.tau_w_max = *std::max_element(profile.tau_w.begin(), profile.tau_w.end());
        profile.tau_w_mean = 0.0;
        for (double val : profile.tau_w) {
            profile.tau_w_mean += val;
        }
        profile.tau_w_mean /= profile.tau_w.size();
        
        return profile;
    }

    /**
     * @brief Validate vortex metrics against reference data
     * 
     * Compares computed metrics with reference values and returns error measures.
     * 
     * @param center_x Computed vortex center X
     * @param center_y Computed vortex center Y
     * @param ref_center_x Reference center X
     * @param ref_center_y Reference center Y
     * @param circulation Computed circulation
     * @param ref_circulation Reference circulation
     * @return VortexValidationResult with error metrics
     */
    static VortexValidationResult validateVortexMetrics(
        double center_x,
        double center_y,
        double ref_center_x,
        double ref_center_y,
        double circulation,
        double ref_circulation)
    {
        VortexValidationResult result;
        
        result.center_x = center_x;
        result.center_y = center_y;
        result.circulation = circulation;
        
        // L2 error in vortex center location
        double dx = center_x - ref_center_x;
        double dy = center_y - ref_center_y;
        result.center_error = std::sqrt(dx * dx + dy * dy);
        
        // Relative circulation error
        if (std::abs(ref_circulation) > 1e-14) {
            result.circulation_error = 100.0 * std::abs(circulation - ref_circulation) / std::abs(ref_circulation);
        }
        
        return result;
    }
};

}
