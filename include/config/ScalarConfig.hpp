#pragma once

#include <string>
#include "../include/linear/Vector.hpp"

namespace config {

// Configuration for generic scalar transport equation
struct ScalarConfig {
    // Scalar identifier
    std::string name = "scalar";  // "temperature", "concentration", etc.
    
    // Transport properties
    linear::Scalar diffusivity = 0.01;    // Γ (thermal diffusivity, mass diffusivity, etc.)
    linear::Scalar turbulent_schmidt = 0.9;  // Sc_t for turbulent diffusion: Γ_t = μ_t / Sc_t
    
    // Physical flags
    bool is_energy = false;       // Special handling for energy (includes dissipation?)
    linear::Scalar heat_capacity = 1.0;   // c_p for energy equation
    linear::Scalar density = 1.0;         // ρ for scalar transport
    
    // Boundary conditions & sources
    bool has_volumetric_source = false;   // Include source term S_φ?
    bool include_dissipation = false;     // For energy: add viscous dissipation Φ?
    bool include_pressure_work = false;   // For energy: add u·∇p term?
    
    // Numerical settings
    int convection_scheme = 0;    // 0=Central, 1=Upwind, 2=QUICK, 3=MUSCL, 4=VanLeer
    int max_scalar_iterations = 50;     // inner loop iterations per time step
    linear::Scalar scalar_tolerance = 1e-6;     // residual tolerance for scalar solver
    
    // Factory methods
    static ScalarConfig temperature() {
        ScalarConfig cfg;
        cfg.name = "temperature";
        cfg.diffusivity = 0.01;     // thermal diffusivity
        cfg.is_energy = true;
        cfg.include_dissipation = true;
        cfg.heat_capacity = 1.0;
        return cfg;
    }
    
    static ScalarConfig concentration() {
        ScalarConfig cfg;
        cfg.name = "concentration";
        cfg.diffusivity = 0.001;    // mass diffusivity
        cfg.is_energy = false;
        cfg.has_volumetric_source = false;
        cfg.turbulent_schmidt = 0.7;   // Often Sc_t < 1 for mass
        return cfg;
    }
    
    static ScalarConfig tracer() {
        ScalarConfig cfg;
        cfg.name = "tracer";
        cfg.diffusivity = 0.0001;   // Very small
        cfg.is_energy = false;
        return cfg;
    }
};

}  // namespace config
