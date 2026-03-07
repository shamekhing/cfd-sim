#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <filesystem>

#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"
#include "../../include/core/UField.hpp"
#include "../../include/core/VField.hpp"
#include "../../include/config/SolverConfig.hpp"
#include "../../include/config/CaseConfigs.hpp"
#include "../../include/solver/StaggeredSIMPLE.hpp"
#include "../../include/solver/ResidualMonitor.hpp"
#include "../../include/physics/TurbulenceKEpsilon.hpp"
#include "../../include/io/VTKWriter.hpp"
#include "../../include/io/CSVWriter.hpp"
#include "../../include/io/PVDWriter.hpp"
#include "../../include/io/VTSReader.hpp"

namespace cases {

/**
 * @brief Turbulent channel flow validation case
 * 
 * Fully developed 2D channel flow driven by constant pressure gradient.
 * Initially laminar, develops turbulence due to instabilities.
 * 
 * Domain: [0, Lx] × [0, 1] (height = 1.0, Ly = 1.0)
 * BCs: No-slip walls at y=0 and y=1
 *      Periodic in x (pressure gradient drives flow)
 * Body force: Constant -dp/dx (in U direction)
 */
inline void runChannelFlow(
    const config::ChannelFlowConfig& cfg,
    const std::string& outputDir = "build/output",
    const std::string& name = "output")
{
    // Domain parameters from config
    const double Lx = cfg.Lx;       // Streamwise extent
    const double Ly = cfg.Ly;       // Channel half-height
    const double nu = cfg.nu();     // Kinematic viscosity
    const double u_bulk = cfg.u_bulk;

    std::cout << "\n=== Turbulent Channel Flow ===" << std::endl;
    std::cout << "  Re_τ (friction Reynolds): " << cfg.Re_tau << std::endl;
    std::cout << "  Domain: [0, " << Lx << "] × [0, " << Ly << "]" << std::endl;
    std::cout << "  Grid: " << cfg.Nx << "x" << cfg.Ny << std::endl;
    std::cout << "  Kinematic viscosity: " << std::scientific << nu << std::defaultfloat << std::endl;

    core::Mesh mesh(cfg.Nx, cfg.Ny, Lx, Ly);

    core::UField u(mesh);
    core::VField v(mesh);
    core::Field p(mesh);

    // Initialize with laminar parabolic profile (Poiseuille flow)
    // u(y) = u_bulk * 6 * y * (1 - y)  [for normalized y in [0,1]]
    // This provides seed for transition to turbulence
    for (int j = 0; j < cfg.Ny; ++j) {
        for (int i = 0; i < cfg.Nx; ++i) {
            double y_norm = mesh.y(j) / Ly;
            u(i, j) = u_bulk * 6.0 * y_norm * (1.0 - y_norm);
            v(i, j) = 0.0;
            p(i, j) = 0.0;
        }
    }

    // Enforce no-slip BC at walls (y=0, y=1)
    for (int i = 0; i < cfg.Nx; ++i) {
        u(i, 0) = 0.0;
        u(i, cfg.Ny - 1) = 0.0;
        v(i, 0) = 0.0;
        v(i, cfg.Ny - 1) = 0.0;
    }

    // Initialize turbulence model with wall functions for high Re
    physics::TurbulenceKEpsilon turbulence;
    config::TurbulenceConfig turb_cfg;
    
    // Use wall functions for high Re channel flow
    if (cfg.useWallFunctions && cfg.Re_tau >= 1000) {
        turb_cfg = config::TurbulenceConfig::withWallFunctions();
        std::cout << "  Wall treatment: logarithmic wall functions" << std::endl;
    } else {
        turb_cfg = config::TurbulenceConfig::standardKEpsilon();
        std::cout << "  Wall treatment: no-slip (low-Re)" << std::endl;
    }
    
    turbulence.initialize(mesh, turb_cfg);

    // Channel-specific configuration from config
    const double alpha_u = cfg.alpha_u;
    const double alpha_v = cfg.alpha_v;
    const double alpha_p = cfg.alpha_p;
    const int maxIter = cfg.maxIter;
    const double tolerance = cfg.tolerance;

    // Pressure gradient to maintain flow (from channel flow theory)
    // τ_w = ρ u_τ² ≈ μ * du/dy|wall; then dp/dx = -2τ_w / h
    const double u_tau = std::sqrt(nu * u_bulk);  // Friction velocity estimate
    const double tau_w = nu * u_tau;               // Wall shear stress
    const double dp_dx = -2.0 * tau_w / Ly;       // Pressure gradient per cell

    std::cout << "  Estimated u_τ: " << std::scientific << u_tau << std::defaultfloat << std::endl;
    std::cout << "  Pressure gradient: " << std::scientific << dp_dx << std::defaultfloat << std::endl;

    const std::string csvOutputDir = outputDir + "/" + name + "/csv";
    const std::string vtsOutputDir = outputDir + "/" + name + "/vts";
    std::filesystem::create_directories(csvOutputDir);
    std::filesystem::create_directories(vtsOutputDir);
    
    // PVD writer for animation with RESUME support
    ::io::PVDWriter pvdWriter(vtsOutputDir, name, cfg.verbosity, outputDir + "/" + name);
    int startIter = pvdWriter.getStartFrame();
    std::string lastFrameFile = pvdWriter.getLastFrameFilename();
    
    // Load field state if resuming
    if (startIter > 0 && !lastFrameFile.empty() && cfg.useAnimation) {
        if (cfg.verbosity >= 1) {
            std::cout << "Resuming from iteration " << startIter << std::endl;
        }
        core::Field u_loaded(mesh);
        core::Field v_loaded(mesh);
        core::Field p_loaded(mesh);
        if (::io::VTSReader::read(lastFrameFile, mesh, u_loaded, v_loaded, p_loaded)) {
            // Copy to staggered fields
            for (int j = 0; j < cfg.Ny; ++j) {
                for (int i = 0; i < cfg.Nx; ++i) {
                    u(i, j) = u_loaded(i, j);
                    v(i, j) = v_loaded(i, j);
                    p(i, j) = p_loaded(i, j);
                }
            }
            solver::StaggeredSIMPLE::applyVelocityBC(u, v, cfg.Nx, cfg.Ny);
        } else if (cfg.verbosity >= 1) {
            std::cout << "Warning: Could not load field state, starting from iter=0\n";
            startIter = 0;
        }
    }

    // SIMPLE iteration loop with turbulence
    auto start_time = std::chrono::high_resolution_clock::now();
    std::ofstream res_file(csvOutputDir + "/residuals.csv");
    res_file << "iter,mass_residual,u_residual,v_residual,k_residual,eps_residual,mass_balance,k_max,eps_max,nu_t_mean,max_yplus,mean_yplus\n";

    core::Field k_old = turbulence.k();
    core::Field eps_old = turbulence.epsilon();

    std::cout << "  Running SIMPLE + k-ε for up to " << maxIter << " iterations..." << std::endl;

    for (int iter = startIter; iter < maxIter; ++iter) {
        // Store old fields for residual computation
        core::UField u_old = u;
        core::VField v_old = v;
        
        // Solve momentum-pressure coupling with body force (pressure gradient)
        // Using SIMPLE solver with dp/dx as source term in U-momentum
        solver::StaggeredSIMPLE::solveSteady(
            mesh, u, v, p,
            nu,
            1,  // Single SIMPLE iteration per outer loop
            alpha_u, alpha_v, alpha_p,
            0,  // No verbosity for inner SIMPLE iteration
            &turbulence,
            dp_dx,  // Source term in U-momentum (pressure gradient)
            0.0);   // No source in V-momentum

        // Compute residuals
        double mass_res = solver::ResidualMonitor::massResidualL2Staggered(mesh, u, v);
        double u_res = (u.raw() - u_old.raw()).norm() / std::max(1e-15, u.raw().norm());
        double v_res = (v.raw() - v_old.raw()).norm() / std::max(1e-15, v.raw().norm());
        double k_res = solver::ResidualMonitor::turbulentKineticEnergyResidualL2(
            turbulence.k(), k_old);
        double eps_res = solver::ResidualMonitor::turbulentDissipationResidualL2(
            turbulence.epsilon(), eps_old);

        k_old = turbulence.k();
        eps_old = turbulence.epsilon();

        // Compute turbulence statistics
        double k_max = 0.0, eps_max = 0.0, nu_t_mean = 0.0;
        for (int idx = 0; idx < turbulence.k().size(); ++idx) {
            k_max = std::max(k_max, turbulence.k()[idx]);
            eps_max = std::max(eps_max, turbulence.epsilon()[idx]);
            nu_t_mean += turbulence.turbulentViscosity()[idx];
        }
        nu_t_mean /= turbulence.k().size();

        // Compute mass balance (∫ ρ u · n dA should be ~0 for periodic BC)
        double mass_balance = 0.0;
        for (int j = 1; j < cfg.Ny - 1; ++j) {
            mass_balance += (u(1, j) - u(cfg.Nx - 2, j));  // Rough check
        }
        mass_balance /= (cfg.Ny - 2);
        
        // Wall y+ diagnostics
        auto [max_yplus, mean_yplus] = solver::ResidualMonitor::wallYPlus(u, v, mesh, nu);

        res_file << iter << "," << mass_res << "," << u_res << "," << v_res << "," 
                 << k_res << "," << eps_res
                 << "," << mass_balance << "," << k_max << "," << eps_max << ","
                 << nu_t_mean << "," << max_yplus << "," << mean_yplus << "\n";

        // Write animation frame every cfg.animationInterval iterations
        // Use iter as "time" so PVD frames have distinct timestep values in ParaView
        if (cfg.useAnimation && iter % cfg.animationInterval == 0) {
            core::Field u_cc = u.cellCenterView(mesh);
            core::Field v_cc = v.cellCenterView(mesh);
            pvdWriter.writeTimestep(iter, static_cast<double>(iter), mesh, u_cc, v_cc, p);
        }

        // Unified progress bar output (consistent format across all cases)
        if ((iter + 1) % 10 == 0 || iter < 3) {
            int progress = (iter * 100) / (maxIter - 1);
            std::cout << "\r  Channel: iter " << iter << "/" << maxIter 
                      << " [" << std::string(progress / 5, '#') 
                      << std::string(20 - progress / 5, '-') << "] " 
                      << progress << "% | mass_res=" << std::scientific << mass_res
                      << std::defaultfloat << std::flush;
        }
        if (iter == maxIter - 1) std::cout << std::endl;

        // Convergence check
        if (mass_res < tolerance && u_res < tolerance && v_res < tolerance && iter > 100) {
            std::cout << "\n  Converged at iteration " << iter << std::endl;
            break;
        }
    }

    res_file.close();

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double>(end_time - start_time).count();
    std::cout << "  Total time: " << std::fixed << std::setprecision(2) << total_time << "s" << std::endl;

    // Extract centerline profiles
    core::Field u_cc = u.cellCenterView(mesh);
    core::Field v_cc = v.cellCenterView(mesh);

    // Compute RMS velocities (turbulent intensities)
    std::vector<double> y_coords, u_mean, u_rms;
    const int probe_i = cfg.Nx / 2;  // x = Lx/2 (middle of domain)

    for (int j = 0; j < cfg.Ny; ++j) {
        y_coords.push_back(mesh.y(j));
        u_mean.push_back(u_cc(probe_i, j));
        // RMS would require time-averaging in true turbulent simulation
        // For now, use turbulent kinetic energy as proxy: u_rms ≈ sqrt(2*k/3)
        double k_val = turbulence.k()[j * cfg.Nx + probe_i];
        u_rms.push_back(std::sqrt(2.0 * k_val / 3.0));
    }

    // Write profiles to CSV
    io::CSVWriter::writeProfile(csvOutputDir + "/u_mean.csv", y_coords, u_mean);
    io::CSVWriter::writeProfile(csvOutputDir + "/u_rms.csv", y_coords, u_rms);

    // Output to VTK
    const std::string vtkFile = outputDir + "/" + name + "/" + name + ".vtk";
    std::filesystem::create_directories(outputDir + "/" + name);
    ::io::VTKWriter::writeWithTurbulence(
        vtkFile,
        mesh,
        u_cc, v_cc, p,
        turbulence.k(),
        turbulence.epsilon(),
        turbulence.turbulentViscosity());

    std::cout << "  VTK output: " << vtkFile << std::endl;
    std::cout << "  Residual history: " << csvOutputDir << "/residuals.csv" << std::endl;
    std::cout << "  Mean velocity profile: " << csvOutputDir << "/u_mean.csv" << std::endl;
    std::cout << "  RMS velocity profile: " << csvOutputDir << "/u_rms.csv" << std::endl;

    // Turbulence statistics summary
    double k_max = 0.0, k_mean = 0.0, eps_max = 0.0;
    for (int idx = 0; idx < turbulence.k().size(); ++idx) {
        k_max = std::max(k_max, turbulence.k()[idx]);
        k_mean += turbulence.k()[idx];
        eps_max = std::max(eps_max, turbulence.epsilon()[idx]);
    }
    k_mean /= turbulence.k().size();

    auto [nu_t_mean, nu_t_max] = solver::ResidualMonitor::turbulentViscosityRatio(
        turbulence.turbulentViscosity(), nu);

    std::cout << "\n  Turbulence Statistics:" << std::endl;
    std::cout << "    k_max: " << std::scientific << k_max << std::defaultfloat << std::endl;
    std::cout << "    k_mean: " << std::scientific << k_mean << std::defaultfloat << std::endl;
    std::cout << "    ε_max: " << std::scientific << eps_max << std::defaultfloat << std::endl;
    std::cout << "    νt/ν mean: " << std::fixed << nu_t_mean << "x" << std::endl;
    std::cout << "    νt/ν max: " << std::fixed << nu_t_max << "x" << std::endl;
    
    // Wall y+ diagnostics
    auto [final_max_yplus, final_mean_yplus] = solver::ResidualMonitor::wallYPlus(u, v, mesh, nu);
    std::cout << "\n  Wall y+ Statistics:" << std::endl;
    std::cout << "    y+_max: " << std::fixed << std::setprecision(2) << final_max_yplus << std::endl;
    std::cout << "    y+_mean: " << std::fixed << std::setprecision(2) << final_mean_yplus << std::endl;
    if (final_max_yplus > 11.25) {
        std::cout << "    (in log-layer region - wall functions appropriate)" << std::endl;
    } else {
        std::cout << "    (in viscous sublayer - low-Re treatment required)" << std::endl;
    }
    
    // Finalize animation
    if (cfg.useAnimation) {
        pvdWriter.finalize();
        pvdWriter.writePVD();
    }
    
    std::cout << "\n=== Channel Flow Complete ===" << std::endl;
}

}  // namespace cases
