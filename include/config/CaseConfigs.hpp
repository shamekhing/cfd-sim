#pragma once

#include <string>
#include "../linear/Vector.hpp"
#include "../io/ColorMap.hpp"

namespace config {

/**
 * @brief Configuration for Lid-Driven Cavity test case
 * 
 * Standard benchmark case with moving lid at top wall
 * Cell-centered test domain: [0,1]x[0,1]
 */
struct LidDrivenCavityConfig {
    // Physical parameters
    linear::Scalar Re = 100.0;
    linear::Scalar Lx = 1.0;
    linear::Scalar Ly = 1.0;
    
    // Grid
    int Nx = 64;
    int Ny = 64;
    
    // Time stepping
    linear::Scalar dt = 1e-2;
    int maxSteadyIter = 500;        // Steady-state iterations
    int numTimeSteps = 1500;          // Transient time steps
    int itersPerTimeStep = 10;
    
    // Relaxation factors (SIMPLE)
    linear::Scalar alpha_u = 0.6;
    linear::Scalar alpha_v = 0.6;
    linear::Scalar alpha_p = 0.3;
    
    // Convergence
    linear::Scalar tolerance = 1e-6;
    int maxLinearIter = 2000;
    linear::Scalar linearTol = 1e-6;
    
    // Animation
    bool useAnimation = true;
    int animationInterval = 1;      // Write every time step
    
    // Validation
    bool validationMode = false;
    
    // Turbulence
    bool enable_turbulence = false;
    
    // Color visualization
    io::ColorScheme colorScheme = io::ColorScheme::BLUE_GREEN_RED;
    bool useColorMap = true;
    
    // Vortex metrics
    bool computeVortexMetrics = true;
    bool outputVorticity = true;
    bool outputStreamfunction = true;
    double vortexDetectionThreshold = 0.05;
    
    int verbosity = 1;
    
    // Compute derived quantities
    linear::Scalar nu() const { return 1.0 / Re; }
};

/**
 * @brief Configuration for High Reynolds Cavity (turbulent regime)
 * 
 * Turbulent cavity flow with k-epsilon model
 * Domain: [0,1]x[0,1]
 */
struct HighReynoldsCavityConfig {
    // Physical parameters
    linear::Scalar Re = 1000.0;
    linear::Scalar Lx = 1.0;
    linear::Scalar Ly = 1.0;
    
    // Grid
    int Nx = 64;
    int Ny = 64;
    
    // Solver iterations
    int maxIter = 800;
    
    // Relaxation factors
    linear::Scalar alpha_u = 0.6;
    linear::Scalar alpha_v = 0.6;
    linear::Scalar alpha_p = 0.3;
    
    // Convergence
    linear::Scalar tolerance = 1e-5;
    int maxLinearIter = 2000;
    linear::Scalar linearTol = 1e-6;
    
    // Time stepping (for potential transient mode)
    linear::Scalar dt = 1e-2;
    int numTimeSteps = 1500;
    
    // Animation
    bool useAnimation = true;
    int animationInterval = 1;      // Write every iteration
    
    // Turbulence model
    bool useWallFunctions = true;   // true if Re >= 1000
    bool enable_turbulence = true;  // Always true for high_re_cavity
    
    int verbosity = 1;
    
    // Compute derived quantities
    linear::Scalar nu() const { return 1.0 / Re; }
};

/**
 * @brief Configuration for Turbulent Channel Flow
 * 
 * Periodic channel with constant pressure gradient
 * Domain: [0, 2π] × [0, 1]
 */
struct ChannelFlowConfig {
    // Physical parameters
    linear::Scalar Re_tau = 5600.0;  // Friction Reynolds number
    linear::Scalar Lx = 2.0 * M_PI;  // Streamwise extent
    linear::Scalar Ly = 1.0;          // Channel half-height
    
    // Grid
    int Nx = 128;
    int Ny = 64;
    
    // Flow initialization
    linear::Scalar u_bulk = 1.0;    // Bulk velocity for Poiseuille profile
    
    // Solver iterations
    int maxIter = 1500;
    
    // Relaxation factors
    linear::Scalar alpha_u = 0.7;
    linear::Scalar alpha_v = 0.7;
    linear::Scalar alpha_p = 0.3;
    
    // Convergence
    linear::Scalar tolerance = 1e-5;
    int maxLinearIter = 2000;
    linear::Scalar linearTol = 1e-6;
    
    // Time stepping
    linear::Scalar dt = 1e-2;
    int numTimeSteps = 100;
    
    // Animation
    bool useAnimation = true;
    int animationInterval = 1;      // Write every iteration
    
    // Turbulence model
    bool useWallFunctions = true;
    
    int verbosity = 1;
    
    // Compute derived quantities
    linear::Scalar nu() const { return 1.0 / Re_tau; }
};

/**
 * @brief Configuration for Energy Transport (MMS validation test)
 * 
 * Tests scalar transport with manufactured solution
 * Domain: [0,1]x[0,1], pure diffusion with exponential decay
 */
struct EnergyTransportConfig {
    // Grid
    int Nx = 64;
    int Ny = 64;
    linear::Scalar Lx = 1.0;
    linear::Scalar Ly = 1.0;
    
    // Time stepping
    linear::Scalar dt = 1e-4;
    int numTimeSteps = 500;
    
    // Thermal diffusivity (scalar diffusion coefficient)
    linear::Scalar alpha_thermal = 0.01;
    
    // Animation
    bool useAnimation = true;
    int animationInterval = 1;      // Write every time step
    
    int verbosity = 1;
};

/**
 * @brief Configuration for Cylinder Flow (vortex shedding dynamics)
 * 
 * Flow around a circular cylinder with Karman vortex shedding
 * Domain: [0, Lx] × [-Ly/2, Ly/2], cylinder at center
 */
struct CylinderFlowConfig {
    // Physical parameters
    linear::Scalar Re = 100.0;      // Reynolds number based on diameter
    linear::Scalar U_inf = 1.0;     // Free-stream velocity
    linear::Scalar radius = 2.0;    // Cylinder radius (diameter D = 2*radius = 4)
    
    // Domain (rectangular, cylinder at center)
    // Sized so that D/dx ≈ 12 on the default grid (adequate for IBM shedding).
    // Blockage D/Ly = 4/24 ≈ 17%, clearance ≈ 2.5D each side.
    linear::Scalar Lx = 40.0;       // Streamwise length  (10D)
    linear::Scalar Ly = 24.0;       // Transverse height  ( 6D)
    
    // Grid – square cells: dx = dy = 1/3
    int Nx = 120;
    int Ny = 72;
    
    // Solver
    int maxIter = 200;              // Per timestep
    linear::Scalar alpha_u = 0.6;
    linear::Scalar alpha_v = 0.6;
    linear::Scalar alpha_p = 0.3;
    
    // Convergence
    linear::Scalar tolerance = 1e-5;
    int maxLinearIter = 2000;
    linear::Scalar linearTol = 1e-6;
    
    // Time stepping (transient for vortex shedding)
    // dt=0.05: CFL ≈ 0.15 on default grid, stable for implicit SIMPLE.
    // numTimeSteps=5000: 250 s physical time ≫ shedding period (~24 s at Re=100).
    linear::Scalar dt = 0.05;
    int numTimeSteps = 5000;        // Need full Strouhal period for validation
    int itersPerTimeStep = 5;       // Inner SIMPLE iters/step (≥5 for deferred-correction convergence)
    
    // Animation
    bool useAnimation = true;
    int animationInterval = 1;      // Write every time step
    
    // Turbulence
    bool enable_turbulence = false; // Use laminar for Re=100, enable for Re>500
    
    int verbosity = 1;
    
    // Compute derived quantities
    linear::Scalar nu() const { return U_inf * (2.0 * radius) / Re; }  // nu = U*D / Re
};

}  // namespace config
