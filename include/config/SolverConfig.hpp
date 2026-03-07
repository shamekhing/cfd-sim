#pragma once

#include "../include/linear/Vector.hpp"
#include "TurbulenceConfig.hpp"

namespace config {

struct SolverConfig {
    // Reynolds number and viscosity
    linear::Scalar Re = 100.0;
    linear::Scalar nu = 0.01;  // 1/Re for Re=100

    // Grid
    int Nx = 64;
    int Ny = 64;
    linear::Scalar Lx = 1.0;
    linear::Scalar Ly = 1.0;

    // Time stepping - GLOBAL DEFAULTS
    linear::Scalar dt = 1e-2;
    int numTimeSteps = 100;  // Global default for transient and steady iterations
    int itersPerTimeStep = 5;

    // Relaxation factors (SIMPLE) - GLOBAL DEFAULTS
    linear::Scalar alpha_u = 0.4;
    linear::Scalar alpha_v = 0.4;
    linear::Scalar alpha_p = 0.2;

    // Linear solver
    int maxLinearIter = 2000;
    linear::Scalar linearTol = 1e-6;

    // Verbosity: 0=quiet, 1=summary, 2=debug
    int verbosity = 1;

    // Animation settings
    std::string animationName = "animation";

    // Steady mode: use SIMPLE::solve for steady-state (meaningful Ghia comparison)
    bool useSteadyMode = true;

    // Validation mode: stricter defaults for benchmark-quality comparisons
    bool validationMode = false;

    // Turbulence modeling
    bool enable_turbulence = false;
    TurbulenceConfig turb_cfg;

    void updateViscosity() {
        nu = 1.0 / Re;
    }

    // Minimal factory - just sets Re and grid, everything else is user-controllable
    static SolverConfig lidDrivenCavity(linear::Scalar re, int nx = 64, int ny = 64) {
        SolverConfig cfg;
        cfg.Re = re;
        cfg.updateViscosity();
        cfg.Nx = nx;
        cfg.Ny = ny;
        return cfg;
    }

    // For animation - user must set useSteadyMode=false, dt, numTimeSteps, etc
    static SolverConfig lidDrivenCavityAnimated(linear::Scalar re, int nx = 64, int ny = 64) {
        SolverConfig cfg = lidDrivenCavity(re, nx, ny);
        cfg.useSteadyMode = false;
        // Everything else uses defaults from struct - user can override
        return cfg;
    }

    static SolverConfig lidDrivenCavityValidation(linear::Scalar re, int nx = 128, int ny = 128) {
        SolverConfig cfg = lidDrivenCavity(re, nx, ny);
        cfg.validationMode = true;
        cfg.useSteadyMode = true;
        cfg.numTimeSteps = 600;
        cfg.alpha_u = 0.6;
        cfg.alpha_v = 0.6;
        cfg.alpha_p = 0.25;
        cfg.verbosity = 1;
        return cfg;
    }
};

}
