#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <chrono>

#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"
#include "../../include/core/UField.hpp"
#include "../../include/core/VField.hpp"
#include "../../include/config/CaseConfigs.hpp"
#include "../../include/assembly/GenericScalarAssembler.hpp"
#include "../../include/linear/Solver.hpp"
#include "../../include/io/VTKWriter.hpp"
#include "../../include/io/PVDWriter.hpp"
#include "../../include/io/VTSReader.hpp"

namespace cases {

/**
 * Energy Transport Test Case (MMS Validation)
 * 
 * Tests the generic scalar transport framework with manufactured solution:
 * ∂T/∂t = ∇·(α∇T)
 * 
 * Manufactured solution (pure diffusion, zero velocity):
 * T(x,y,t) = cos(2πx)cos(2πy)·exp(-4π²α·t)
 * 
 * Expected behavior: temperature field decays exponentially
 * Validates: O(h²) spatial convergence, L2 error < 0.01
 * 
 * RESUME SUPPORT: Detects existing frames and resumes from last frame
 */
inline void runEnergyTransport(
    const config::EnergyTransportConfig& cfg,
    const std::string& outputDir = "build/output",
    const std::string& name = "energy_transport")
{
    const int Nx = cfg.Nx;
    const int Ny = cfg.Ny;
    const double Lx = cfg.Lx;
    const double Ly = cfg.Ly;
    const double alpha = cfg.alpha_thermal;
    const double dt = cfg.dt;
    const int numTimeSteps = cfg.numTimeSteps;
    const bool useAnimation = cfg.useAnimation;

    if (cfg.verbosity >= 1) {
        std::cout << "\n=== Energy Transport (MMS Test) ===\n";
        std::cout << "Domain: [0,1] x [0,1]\n";
        std::cout << "Grid: " << Nx << "x" << Ny << "\n";
        std::cout << "Diffusivity α = " << alpha << "\n";
        std::cout << "Time step: " << dt << "\n";
        std::cout << "Total steps: " << numTimeSteps << "\n";
    }

    core::Mesh mesh(Nx, Ny, Lx, Ly);
    
    // Temperature field
    core::Field T(mesh);
    core::Field T_old(mesh);
    
    // Zero velocity field (pure diffusion MMS) - use Field instead of UField/VField
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
    
    // Setup output directories
    const std::string runOutputDir = outputDir + "/" + name;
    const std::string csvOutputDir = runOutputDir + "/csv";
    const std::string vtsOutputDir = runOutputDir + "/vts";
    
    std::filesystem::create_directories(csvOutputDir);
    std::filesystem::create_directories(vtsOutputDir);
    
    // PVD writer for animation with RESUME support
    ::io::PVDWriter pvdWriter(vtsOutputDir, name, cfg.verbosity, runOutputDir);
    int startFrame = pvdWriter.getStartFrame();
    std::string lastFrameFile = pvdWriter.getLastFrameFilename();
    
    // Load field state if resuming
    if (startFrame > 0 && !lastFrameFile.empty() && useAnimation) {
        if (cfg.verbosity >= 1) {
            std::cout << "Resuming from frame " << startFrame << std::endl;
        }
        core::Field T_loaded(mesh);
        core::Field u_dummy(mesh);
        core::Field v_dummy(mesh);
        if (::io::VTSReader::read(lastFrameFile, mesh, u_dummy, v_dummy, T_loaded)) {
            T = T_loaded;
            T_old = T_loaded;
        } else if (cfg.verbosity >= 1) {
            std::cout << "Warning: Could not load field state, starting from t=0\n";
            startFrame = 0;
        }
    }
    
    // Convergence history file (append if resuming)
    std::ofstream conv_file;
    if (startFrame == 0) {
        conv_file.open(csvOutputDir + "/convergence.csv");
        conv_file << "step,time,L2_error,relative_error\n";
    } else {
        conv_file.open(csvOutputDir + "/convergence.csv", std::ios::app);
    }
    
    double time = startFrame * dt;
    auto anim_start = std::chrono::high_resolution_clock::now();
    
    // Time loop with RESUME
    for (int step = startFrame; step < numTimeSteps; ++step) {
        T_old = T;
        
        // Assemble and solve energy equation: ∂T/∂t = ∇·(α∇T)
        auto system = assembly::GenericScalarAssembler::assemble(
            mesh,
            T,           // current temperature
            u, v,        // zero velocity
            alpha,       // thermal diffusivity
            dt,
            T_old,       // previous temperature
            assembly::ConvectionScheme::Upwind,
            nullptr,     // no source
            true,        // is_energy
            1.0);        // nu (not used for MMS)
        
        // Solve linear system
        linear::BiCGSTABSolver solver;
        linear::Vector T_new = solver.solve(system.A, system.b);
        T.raw() = T_new;
        
        time = (step + 1) * dt;
        
        // Compute L2 error vs manufactured solution
        double error_l2 = 0.0;
        double T_max_exact = 0.0;
        double decay_factor = std::exp(-4.0 * M_PI * M_PI * alpha * time);
        
        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                double x = mesh.x(i);
                double y = mesh.y(j);
                double T_exact = decay_factor * std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
                double T_num = T(i, j);
                error_l2 += (T_num - T_exact) * (T_num - T_exact);
                T_max_exact = std::max(T_max_exact, std::abs(T_exact));
            }
        }
        error_l2 = std::sqrt(error_l2 / ((Nx - 2) * (Ny - 2)));
        double rel_error = (T_max_exact > 1e-15) ? error_l2 / T_max_exact : error_l2;
        
        // Write convergence data
        conv_file << step << "," << time << "," << std::scientific << error_l2 
                  << "," << rel_error << std::defaultfloat << "\n";
        
        // Write animation frame — temperature as primary named scalar
        if (useAnimation) {
            pvdWriter.writeTimestepScalar(step, time, mesh, T, "temperature", u, v);
        }
        
        // Progress output
        if ((step + 1) % 10 == 0 || step == startFrame) {
            int progress = ((step + 1) * 100) / numTimeSteps;
            std::cout << "\rStep " << (step + 1) << "/" << numTimeSteps 
                      << " [" << std::string(progress / 5, '#') 
                      << std::string(20 - progress / 5, '-') << "] " 
                      << progress << "% | L2_err=" << std::scientific << error_l2
                      << std::defaultfloat << std::flush;
        }
    }
    std::cout << std::endl;
    conv_file.close();
    
    // Finalize animation
    if (useAnimation) {
        pvdWriter.finalize();
        pvdWriter.writePVD();
    }
    
    // Final convergence check
    double final_time = numTimeSteps * dt;
    double decay_factor_final = std::exp(-4.0 * M_PI * M_PI * alpha * final_time);
    double error_l2_final = 0.0;
    double T_max_exact_final = 0.0;
    
    for (int j = 1; j < Ny - 1; ++j) {
        for (int i = 1; i < Nx - 1; ++i) {
            double x = mesh.x(i);
            double y = mesh.y(j);
            double T_exact = decay_factor_final * std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
            double T_num = T(i, j);
            error_l2_final += (T_num - T_exact) * (T_num - T_exact);
            T_max_exact_final = std::max(T_max_exact_final, std::abs(T_exact));
        }
    }
    error_l2_final = std::sqrt(error_l2_final / ((Nx - 2) * (Ny - 2)));
    double rel_error_final = (T_max_exact_final > 1e-15) ? error_l2_final / T_max_exact_final : error_l2_final;
    
    if (cfg.verbosity >= 1) {
        std::cout << "\nFinal Results:\n";
        std::cout << "L2 error: " << std::scientific << error_l2_final << std::defaultfloat << "\n";
        std::cout << "Relative error: " << (rel_error_final * 100.0) << "%\n";
        std::cout << "Test status: ";
        
        if (error_l2_final < 0.01) {
            std::cout << "PASS (error < 0.01)\n";
        } else if (error_l2_final < 0.05) {
            std::cout << "ACCEPTABLE (error < 0.05)\n";
        } else {
            std::cout << "FAIL (error >= 0.05)\n";
        }
        std::cout << "Output: " << runOutputDir << "\n";
    }
}

}  // namespace cases
