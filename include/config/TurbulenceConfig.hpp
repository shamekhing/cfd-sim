#pragma once

#include <string>
#include <cmath>
#include "../include/linear/Vector.hpp"

namespace config {

/**
 * @brief Configuration for k-epsilon turbulence model
 * 
 * Standard k-epsilon model (Launder & Spalding, 1974):
 *   dk/dt + u·∇k = ∇·[(ν + νt/σk)∇k] + Pk - ρε
 *   dε/dt + u·∇ε = ∇·[(ν + νt/σε)∇ε] + C1ε(ε/k)Pk - C2ε ρ(ε²/k)
 *   νt = Cμ k²/ε
 */
struct TurbulenceConfig
{
    // Model constants (standard k-epsilon)
    linear::Scalar C_mu = 0.09;          // Turbulent viscosity coefficient
    linear::Scalar sigma_k = 1.00;       // Prandtl number for k
    linear::Scalar sigma_epsilon = 1.30; // Prandtl number for epsilon
    linear::Scalar C_1epsilon = 1.44;    // Production coefficient in epsilon equation
    linear::Scalar C_2epsilon = 1.92;    // Dissipation coefficient in epsilon equation
    
    // Numerical parameters
    int inner_iterations = 2;    // Inner k-epsilon iterations per SIMPLE step
    linear::Scalar relax_k = 0.7;        // Under-relaxation for k
    linear::Scalar relax_epsilon = 0.7;  // Under-relaxation for epsilon
    linear::Scalar min_k = 1e-10;        // Minimum k to prevent division by zero
    linear::Scalar min_epsilon = 1e-10;  // Minimum epsilon to prevent division by zero
    
    // Boundary condition parameters
    enum WallBCType { NoSlip, WallFunction, Symmetry };
    WallBCType wall_bc_type = NoSlip;
    
    // Inlet turbulence specification
    linear::Scalar turbulence_intensity = 0.05;  // I = u'/U (typical: 1-5%)
    linear::Scalar turbulent_length_scale = 0.1; // Lt (estimated from geometry)
    
    /**
     * @brief Factory: Standard k-epsilon for lid-driven cavity
     */
    static TurbulenceConfig standardKEpsilon()
    {
        TurbulenceConfig cfg;
        cfg.wall_bc_type = NoSlip;
        cfg.inner_iterations = 2;
        return cfg;
    }
    
    /**
     * @brief Factory: k-epsilon with wall functions (for high Re)
     */
    static TurbulenceConfig withWallFunctions()
    {
        TurbulenceConfig cfg;
        cfg.wall_bc_type = WallFunction;
        cfg.inner_iterations = 2;
        return cfg;
    }
    
    /**
     * @brief Compute inlet k from turbulence intensity and reference velocity
     * k_inlet = 1.5 * (I * U_ref)²
     */
    linear::Scalar computeInletK(linear::Scalar U_ref) const
    {
        return 1.5 * (turbulence_intensity * U_ref) * (turbulence_intensity * U_ref);
    }
    
    /**
     * @brief Compute inlet epsilon from k and turbulent length scale
     * ε_inlet = C_mu^(3/4) * k^(3/2) / L_t
     */
    linear::Scalar computeInletEpsilon(linear::Scalar k_inlet) const
    {
        return std::pow(C_mu, 0.75) * std::pow(k_inlet, 1.5) / turbulent_length_scale;
    }
};

} // namespace config
