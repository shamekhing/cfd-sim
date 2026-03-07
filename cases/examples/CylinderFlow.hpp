#pragma once

#include <iostream>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>

#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"
#include "../../include/core/UField.hpp"
#include "../../include/core/VField.hpp"
#include "../../include/core/Indexing.hpp"
#include "../../include/config/SolverConfig.hpp"
#include "../../include/config/CaseConfigs.hpp"
#include "../../include/solver/StaggeredSIMPLE.hpp"
#include "../../include/solver/ResidualMonitor.hpp"
#include "../../include/io/VTKWriter.hpp"
#include "../../include/io/CSVWriter.hpp"
#include "../../include/io/PVDWriter.hpp"
#include "../../include/io/VTSReader.hpp"

namespace cases {

/**
 * @brief Cylinder Flow Case (Vortex Shedding Validation)
 * 
 * Flow around a circular cylinder for testing vortex dynamics and Strouhal number
 * Domain: rectangular [0, Lx] × [-Ly/2, Ly/2]
 * Cylinder: centered at (Lc, 0) with radius r
 * 
 * Boundary conditions:
 * - Inlet (left, x=0): u = U_inf, v = 0 (uniform)
 * - Outlet (right, x=Lx): ∂u/∂x = 0, v = 0 (zero-gradient)
 * - Top/bottom walls: y = ±Ly/2, slip (free-stream)
 * - Cylinder: r from center = exact circle, no-slip
 * 
 * Physics:
 * - Laminar vortex shedding at Re=100: f_shedding ≈ U/D * St with St ≈ 0.164
 * - Strouhal number: St = f·D / U where f is shedding frequency, D is diameter
 * 
 * Output:
 * - VTS files for visualization (velocity magnitude, pressure)
 * - CSV with timestamp, lift coefficient, drag coefficient
 * - Residual history
 */
inline void runCylinderFlow(
    const config::CylinderFlowConfig& cfg,
    const std::string& outputDir = "build/output",
    const std::string& name = "cylinder_flow")
{
    const int Nx = cfg.Nx;
    const int Ny = cfg.Ny;
    const double Lx = cfg.Lx;
    const double Ly = cfg.Ly;
    const double Re = cfg.Re;
    const double U_inf = cfg.U_inf;
    const double r = cfg.radius;
    const double nu = cfg.nu();
    const double dt = cfg.dt;
    const int numTimeSteps = cfg.numTimeSteps;
    
    // Cylinder center location (roughly 1/4 of domain from inlet)
    const double cx = Lx / 4.0;
    const double cy = 0.0;
    
    const std::string runOutputDir = outputDir + "/" + name;
    const std::string vtsOutputDir = runOutputDir + "/vts";
    const std::string csvOutputDir = runOutputDir + "/csv";
    std::filesystem::create_directories(vtsOutputDir);
    std::filesystem::create_directories(csvOutputDir);
    
    if (cfg.verbosity >= 1) {
        std::cout << "\n=== Cylinder Flow (Vortex Shedding) ===\n";
        std::cout << "Domain: [0, " << Lx << "] × [" << -Ly/2.0 << ", " << Ly/2.0 << "]\n";
        std::cout << "Cylinder: center (" << cx << ", " << cy << "), radius " << r << "\n";
        std::cout << "Reynolds number Re = " << Re << ", ν = " << nu << "\n";
        std::cout << "Free-stream velocity U_∞ = " << U_inf << "\n";
        std::cout << "Grid: " << Nx << " × " << Ny << "\n";
        std::cout << "Time step: " << dt << ", Total steps: " << numTimeSteps << "\n";
        std::cout << "Expected Strouhal (Re=100): St ≈ 0.164\n";
        std::cout << "Expected shedding frequency: f ≈ " << (U_inf / (2.0 * r) * 0.164) << " Hz\n";
        std::cout << "Period: T ≈ " << (2.0 * r / U_inf / 0.164) << " s\n\n";
    }
    
    // Create mesh (note: this creates a regular Cartesian grid; cylinder is handled via BCs)
    core::Mesh mesh(Nx, Ny, Lx, Ly);
    
    // Initialize fields
    core::UField u(mesh);
    core::VField v(mesh);
    core::Field p(mesh);
    
    u.setZero();
    v.setZero();
    p.setZero();
    
    // Setup PVD writer for animation with RESUME support
    ::io::PVDWriter pvdWriter(vtsOutputDir, name, cfg.verbosity, runOutputDir);
    int startFrame = pvdWriter.getStartFrame();
    std::string lastFrameFile = pvdWriter.getLastFrameFilename();
    
    // Load field state if resuming
    if (startFrame > 0 && !lastFrameFile.empty() && cfg.useAnimation) {
        if (cfg.verbosity >= 1) {
            std::cout << "Resuming from frame " << startFrame << std::endl;
        }
        core::Field u_loaded(mesh);
        core::Field v_loaded(mesh);
        core::Field p_loaded(mesh);
        if (::io::VTSReader::read(lastFrameFile, mesh, u_loaded, v_loaded, p_loaded)) {
            // Load into staggered fields (simplified: direct copy)
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    u(i, j) = u_loaded(i, j);
                    v(i, j) = v_loaded(i, j);
                    p(i, j) = p_loaded(i, j);
                }
            }
        } else if (cfg.verbosity >= 1) {
            std::cout << "Warning: Could not load field state, starting from t=0\n";
            startFrame = 0;
        }
    }
    
    // Initial condition: uniform free-stream everywhere (skip if resuming)
    if (startFrame == 0) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = mesh.x(i);
                double y = mesh.y(j);
                
                // Check if point is inside or very close to cylinder
                double dx = x - cx;
                double dy = y - cy;
                double dist = std::sqrt(dx*dx + dy*dy);
                
                if (dist > r + mesh.dx()) {
                    u(i, j) = U_inf;
                }
            }
        }
    }
    
    // Solver configuration
    config::SolverConfig solverCfg;
    solverCfg.nu = nu;
    solverCfg.dt = dt;
    solverCfg.alpha_u = cfg.alpha_u;
    solverCfg.alpha_v = cfg.alpha_v;
    solverCfg.alpha_p = cfg.alpha_p;
    solverCfg.verbosity = cfg.verbosity;
    
    // Output files
    std::ofstream res_file(csvOutputDir + "/residuals.csv");
    res_file << "step,time,mass_residual,u_residual,v_residual,p_residual\n";
    
    std::ofstream force_file(csvOutputDir + "/forces.csv");
    force_file << "step,time,lift_coeff,drag_coeff\n";
    
    double time = startFrame * dt;
    auto anim_start = std::chrono::high_resolution_clock::now();

    // Time stepping loop with RESUME support
    for (int step = startFrame; step < numTimeSteps; ++step) {
        // Apply boundary conditions
        // Inlet (left boundary, x=0): u = U_inf
        for (int j = 0; j < Ny; ++j) {
            u(0, j) = U_inf;
        }
        
        // Apply no-slip on cylinder surface
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = mesh.x(i) + mesh.dx() / 2.0;
                double y = mesh.y(j);
                
                double dx = x - cx;
                double dy = y - cy;
                double dist = std::sqrt(dx*dx + dy*dy);
                
                // Apply no-slip condition within cylinder region
                if (dist <= r) {
                    u(i, j) = 0.0;
                }
                
                // V-velocity at (i, j+0.5)
                double yv = mesh.y(j) + mesh.dy() / 2.0;
                double dyv = yv - cy;
                double distv = std::sqrt(dx*dx + dyv*dyv);
                if (distv <= r) {
                    v(i, j) = 0.0;
                }
            }
        }
        
        // Solve one SIMPLE iteration for steady state or pseudo-transient
        double mass_res = 0.0, u_res = 0.0, v_res = 0.0, p_res = 0.0;
        solver::StaggeredSIMPLE::solveSteady(
            mesh, u, v, p, nu,
            20,  // iterations per timestep
            solverCfg.alpha_u, solverCfg.alpha_v, solverCfg.alpha_p,
            0);  // suppress per-iteration output; progress bar handles display
        
        // Re-apply cylinder boundary conditions after solve
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = mesh.x(i) + mesh.dx() / 2.0;
                double y = mesh.y(j);
                
                double dx = x - cx;
                double dy = y - cy;
                double dist = std::sqrt(dx*dx + dy*dy);
                
                if (dist <= r) {
                    u(i, j) = 0.0;
                }
                
                double yv = mesh.y(j) + mesh.dy() / 2.0;
                double dyv = yv - cy;
                double distv = std::sqrt(dx*dx + dyv*dyv);
                if (distv <= r) {
                    v(i, j) = 0.0;
                }
            }
        }
        
        time += dt;
        
        // Simplified force calculation (pressure integral around cylinder)
        double lift_coeff = 0.0;
        double drag_coeff = 0.0;
        
        // Very rough estimate: collect pressure samples on cylinder
        int sample_count = 0;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = mesh.x(i);
                double y = mesh.y(j);
                
                double dx = x - cx;
                double dy = y - cy;
                double dist = std::sqrt(dx*dx + dy*dy);
                
                // Sample pressure near cylinder surface
                if (dist > r && dist < r + mesh.dx() * 2.0) {
                    double normal_x = dx / dist;
                    double normal_y = dy / dist;
                    
                    drag_coeff += p(i, j) * normal_x;
                    lift_coeff += p(i, j) * normal_y;
                    sample_count++;
                }
            }
        }
        
        if (sample_count > 0) {
            double diam = 2.0 * r;
            double ref_area = diam * 1.0;  // Assume unit depth
            double dyn_pressure = 0.5 * U_inf * U_inf;
            
            drag_coeff /= (sample_count * dyn_pressure * ref_area);
            lift_coeff /= (sample_count * dyn_pressure * ref_area);
        }
        
        {
            auto current_time = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(current_time - anim_start).count();
            double avg_time = elapsed / (step - startFrame + 1);
            int remaining = numTimeSteps - step - 1;
            double eta = avg_time * remaining;
            int progress = ((step + 1) * 100) / numTimeSteps;
            std::string prefix = (startFrame > 0) ? "  Resuming: " : "  Animation: ";
            std::cout << "\r" << prefix << (step + 1) << "/" << numTimeSteps
                      << " [" << std::string(progress / 5, '#')
                      << std::string(20 - progress / 5, '-') << "] "
                      << progress << "% | "
                      << std::fixed << std::setprecision(1)
                      << "T:" << elapsed << "s | ETA:" << eta << "s"
                      << " | drag=" << std::setprecision(3) << drag_coeff
                      << " lift=" << lift_coeff
                      << std::flush;
        }
        
        res_file << step + 1 << "," << time << "," << mass_res << "," << u_res << "," << v_res << "," << p_res << "\n";
        force_file << step + 1 << "," << time << "," << lift_coeff << "," << drag_coeff << "\n";
        
        // Write VTS for animation using PVDWriter (with cell-center view)
        if (cfg.useAnimation) {
            core::Field u_cc = u.cellCenterView(mesh);
            core::Field v_cc = v.cellCenterView(mesh);
            pvdWriter.writeTimestep(step, time, mesh, u_cc, v_cc, p);
        }
    }
    
    std::cout << std::endl;

    res_file.close();
    force_file.close();
    
    // Finalize animation
    if (cfg.useAnimation) {
        pvdWriter.finalize();
        pvdWriter.writePVD();
    }
    
    if (cfg.verbosity >= 1) {
        std::cout << "\n=== Cylinder Flow Complete ===\n";
        std::cout << "Output written to: " << vtsOutputDir << "/\n";
        std::cout << "Forces: " << csvOutputDir << "/forces.csv\n";
        std::cout << "Residuals: " << csvOutputDir << "/residuals.csv\n";
    }
}

}
