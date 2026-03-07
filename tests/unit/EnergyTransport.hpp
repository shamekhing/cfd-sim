#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "../include/core/UField.hpp"
#include "../include/core/VField.hpp"
#include "../include/config/SolverConfig.hpp"
#include "../include/config/ScalarConfig.hpp"
#include "../include/physics/ScalarTransport.hpp"
#include "../include/assembly/GenericScalarAssembler.hpp"
#include "../include/linear/Solver.hpp"
#include "../include/io/Logger.hpp"

namespace cases {

/**
 * Energy Transport Test Case
 * 
 * Tests the generic scalar transport framework with energy equation:
 * ∂T/∂t + u·∇T = ∇·(α∇T) + Φ
 * 
 * Where Φ = 2μ(S·S) is viscous dissipation
 * 
 * Manufactured solution:
 * T(x,y,t) = cos(2πx)cos(2πy)·exp(-4π²α·t)
 * 
 * Expected behavior: temperature field decays exponentially over time
 */
inline void runEnergyTransport(
    const std::string& outputDir = "build/output",
    const std::string& name = "energy_transport")
{
    const int Nx = 64;
    const int Ny = 64;
    const double Lx = 1.0;
    const double Ly = 1.0;

    core::Mesh mesh(Nx, Ny, Lx, Ly);
    
    // Energy field (temperature)
    core::Field T(mesh);
    core::Field T_old(mesh);
    
    // Thermal diffusivity (simplified: no velocity)
    double alpha = 0.01;
    
    // Velocity field (zero for simplified test, or use manufactured field)
    core::Field u(mesh);
    core::Field v(mesh);
    u.setZero();
    v.setZero();
    
    // Initialize with manufactured solution at t=0
    // T(x,y,0) = cos(2πx)cos(2πy)
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double x = mesh.x(i);
            double y = mesh.y(j);
            T(i, j) = std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
        }
    }
    T_old = T;
    
    // Scalar configuration for energy
    config::ScalarConfig cfg = config::ScalarConfig::temperature();
    cfg.diffusivity = alpha;
    cfg.include_dissipation = false;  // No dissipation in this manufactured solution
    cfg.convection_scheme = 1;  // Upwind for stability
    
    // Time stepping
    double dt = 1e-4;
    int num_steps = 100;
    double t = 0.0;
    
    // Exact solution at final time
    double t_final = t + num_steps * dt;
    double decay_factor = std::exp(-4.0 * M_PI * M_PI * alpha * t_final);
    
    std::cout << "\n=== Energy Transport Test Case ===\n";
    std::cout << "Domain: [0,1] x [0,1]\n";
    std::cout << "Grid: " << Nx << " x " << Ny << "\n";
    std::cout << "Thermal diffusivity α = " << alpha << "\n";
    std::cout << "Time step: " << dt << "\n";
    std::cout << "Number of steps: " << num_steps << "\n";
    std::cout << "Total time: " << t_final << "\n";
    std::cout << "Expected decay factor: " << decay_factor << "\n\n";
    
    // Time loop
    for (int step = 0; step < num_steps; ++step) {
        T_old = T;
        
        // Assemble and solve energy equation: ∂T/∂t = ∇·(α∇T)
        auto system = assembly::GenericScalarAssembler::assemble(
            mesh,
            T,       // current temperature
            u, v,    // zero velocity
            alpha,   // thermal diffusivity
            dt,
            T_old,   // previous temperature
            assembly::ConvectionScheme::Upwind,
            nullptr, // no source
            true,    // is_energy (but no dissipation for MMS)
            0.01);   // nu (not used if no dissipation)
        
        // Solve
        linear::Vector T_new = linear::BiCGSTABSolver().solve(system.A, system.b);
        T.raw() = T_new;
        
        t += dt;
        
        if ((step + 1) % 20 == 0) {
            std::cout << "Step " << (step + 1) << " / " << num_steps
                      << ", t = " << std::setprecision(6) << t << "\n";
        }
    }
    
    // Compute L2 error vs manufactured solution
    double error_l2 = 0.0;
    double T_max_exact = 0.0;
    
    for (int j = 1; j < Ny - 1; ++j) {
        for (int i = 1; i < Nx - 1; ++i) {
            double x = mesh.x(i);
            double y = mesh.y(j);
            
            // Exact solution at final time
            double T_exact = decay_factor * std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
            double T_num = T(i, j);
            
            error_l2 += (T_num - T_exact) * (T_num - T_exact);
            T_max_exact = std::max(T_max_exact, std::abs(T_exact));
        }
    }
    error_l2 = std::sqrt(error_l2 / ((Nx - 2) * (Ny - 2)));
    
    std::cout << "\n=== Results ===\n";
    std::cout << "L2 error: " << std::setprecision(8) << error_l2 << "\n";
    std::cout << "Relative error: " << (error_l2 / T_max_exact) * 100.0 << " %\n";
    std::cout << "Test status: ";
    
    if (error_l2 < 0.01) {
        std::cout << "PASS (error < 0.01)\n";
    } else if (error_l2 < 0.05) {
        std::cout << "ACCEPTABLE (error < 0.05)\n";
    } else {
        std::cout << "FAIL (error >= 0.05)\n";
    }
    
    // Write field for visualization
    const std::string csvFile = outputDir + "/" + name + "/" + name + ".csv";
    std::filesystem::create_directories(outputDir + "/" + name);
    std::ofstream out(csvFile);
    out << "x,y,T_numerical,T_exact\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double x = mesh.x(i);
            double y = mesh.y(j);
            double T_exact = decay_factor * std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
            out << x << "," << y << "," << T(i, j) << "," << T_exact << "\n";
        }
    }
    out.close();
    
    std::cout << "Output written to: " << csvFile << "\n";
}

}  // namespace cases
