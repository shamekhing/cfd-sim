#pragma once

#include <iostream>
#include <cmath>
#include <filesystem>
#include <fstream>

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
#include "../../tests/validation/GhiaReference.hpp"

namespace cases {

/**
 * @brief High Reynolds number lid-driven cavity with turbulence
 * 
 * Test case for validating k-epsilon model at Re=1000 (turbulent regime)
 * Compares results with literature and lower Reynolds numbers.
 */
inline void runHighReynoldsCavity(
    const config::HighReynoldsCavityConfig& cfg,
    const std::string& outputDir = "build/output",
    const std::string& name = "output")
{
    if (cfg.Re < 500) {
        std::cerr << "HighReynoldsCavity: Re should be >= 500 for turbulent effects. Got Re=" << cfg.Re << std::endl;
        return;
    }
    
    const double nu = cfg.nu();
    std::cout << "\n=== High Reynolds Cavity (Turbulent) ===" << std::endl;
    std::cout << "  Re = " << cfg.Re << ", Grid: " << cfg.Nx << "x" << cfg.Ny << std::endl;
    std::cout << "  Turbulence model: k-epsilon" << std::endl;

    core::Mesh mesh(cfg.Nx, cfg.Ny, cfg.Lx, cfg.Ly);

    core::UField u(mesh);
    core::VField v(mesh);
    core::Field p(mesh);

    u.setZero();
    v.setZero();
    p.setZero();
    solver::StaggeredSIMPLE::applyVelocityBC(u, v, cfg.Nx, cfg.Ny);

    // Initialize turbulence model
    physics::TurbulenceKEpsilon turbulence;
    config::TurbulenceConfig turb_cfg;
    
    // Use wall functions for high Re where y+ is expected to be > 11.25
    if (cfg.useWallFunctions && cfg.Re >= 1000) {
        turb_cfg = config::TurbulenceConfig::withWallFunctions();
        std::cout << "  Wall treatment: logarithmic wall functions" << std::endl;
    } else {
        turb_cfg = config::TurbulenceConfig::standardKEpsilon();
        std::cout << "  Wall treatment: no-slip (low-Re)" << std::endl;
    }
    
    turbulence.initialize(mesh, turb_cfg);

    // Create output directory
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

    // Residual history files
    std::ofstream res_file(csvOutputDir + "/residuals.csv");
    res_file << "iter,u_residual,v_residual,p_residual,k_residual,eps_residual,mass_residual,max_yplus,mean_yplus\n";

    // Convergence settings from config
    const int maxIter = cfg.maxIter;
    const linear::Scalar alpha_u = cfg.alpha_u;
    const linear::Scalar alpha_v = 0.5;
    const linear::Scalar alpha_p = 0.2;

    std::cout << "  Relaxation: α_u=" << alpha_u << ", α_v=" << alpha_v << ", α_p=" << alpha_p << std::endl;
    std::cout << "  Running SIMPLE + k-ε for up to " << maxIter << " iterations..." << std::endl;

    // Manual SIMPLE iteration with residual tracking
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int iter = startIter; iter < maxIter; ++iter) {
        core::UField u_old = u;
        core::VField v_old = v;

        // Solve momentum and pressure with unified progress bar
        solver::StaggeredSIMPLE::solveSteady(mesh, u, v, p, nu, 1,
            alpha_u, alpha_v, alpha_p, 0, &turbulence);  // Verbosity=0, outer loop shows progress

        // Compute residuals every 10 iterations
        if (iter % 10 == 0) {
            double mass_res = solver::ResidualMonitor::massResidualL2Staggered(mesh, u, v);
            
            // Field residuals
            double u_res = (u.raw() - u_old.raw()).norm() / std::max(1e-15, u.raw().norm());
            double v_res = (v.raw() - v_old.raw()).norm() / std::max(1e-15, v.raw().norm());
            double p_res = 0.0;  // Placeholder: proper p residual tracking in Phase 2.3

            // Turbulence residuals from k-ε subsystem
            double k_res = turbulence.lastKResidual();
            double eps_res = turbulence.lastEpsilonResidual();
            
            // Wall y+ diagnostics
            auto [max_yplus, mean_yplus] = solver::ResidualMonitor::wallYPlus(u, v, mesh, nu);

            res_file << iter << "," << u_res << "," << v_res << "," << p_res << ","
                     << k_res << "," << eps_res << "," << mass_res << ","
                     << max_yplus << "," << mean_yplus << "\n";

            // Write animation frame every cfg.animationInterval iterations
            // Use iter as "time" so PVD frames have distinct timestep values in ParaView
            if (cfg.useAnimation && iter % cfg.animationInterval == 0) {
                core::Field u_cc = u.cellCenterView(mesh);
                core::Field v_cc = v.cellCenterView(mesh);
                pvdWriter.writeTimestep(iter, static_cast<double>(iter), mesh, u_cc, v_cc, p);
            }

            // Unified progress bar output (consistent format across all cases)
            int progress = (iter * 100) / (maxIter - 1);
            std::cout << "\r  HighRe Cavity: iter " << iter << "/" << maxIter 
                      << " [" << std::string(progress / 5, '#') 
                      << std::string(20 - progress / 5, '-') << "] " 
                      << progress << "% | mass_res=" << std::scientific << mass_res
                      << std::defaultfloat << std::flush;
            if (iter == maxIter - 1) std::cout << std::endl;

            // Convergence check
            if (mass_res < 1e-4 && iter > 50) {
                std::cout << "\n  Converged at iteration " << iter << " (mass residual < 1e-4)" << std::endl;
                break;
            }
        }
    }
    res_file.close();
    std::cout << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double>(end_time - start_time).count();
    std::cout << "  Total time: " << std::fixed << std::setprecision(2) << total_time << "s" << std::endl;

    // Output results
    double mass_res = solver::ResidualMonitor::massResidualL2Staggered(mesh, u, v);
    std::cout << "  Mass residual L2: " << std::scientific << mass_res << std::defaultfloat << std::endl;

    core::Field u_cc = u.cellCenterView(mesh);
    core::Field v_cc = v.cellCenterView(mesh);

    double maxVelInterior = 0.0;
    for (int j = 1; j < cfg.Ny - 1; ++j) {
        for (int i = 1; i < cfg.Nx - 1; ++i) {
            double mag = std::sqrt(u_cc(i, j) * u_cc(i, j) + v_cc(i, j) * v_cc(i, j));
            if (mag > maxVelInterior) maxVelInterior = mag;
        }
    }
    std::cout << "  Max velocity magnitude (interior): " << maxVelInterior << std::endl;

    // Turbulence statistics
    const auto& k_field = turbulence.k();
    const auto& eps_field = turbulence.epsilon();
    const auto& nu_t_field = turbulence.turbulentViscosity();

    double k_max = 0.0, k_mean = 0.0;
    double eps_max = 0.0;
    double nu_t_max = 0.0, nu_t_mean = 0.0;

    for (int idx = 0; idx < k_field.size(); ++idx) {
        k_max = std::max(k_max, k_field[idx]);
        k_mean += k_field[idx];
        eps_max = std::max(eps_max, eps_field[idx]);
        nu_t_max = std::max(nu_t_max, nu_t_field[idx]);
        nu_t_mean += nu_t_field[idx];
    }
    k_mean /= k_field.size();
    nu_t_mean /= nu_t_field.size();

    std::cout << "  Turbulence statistics:" << std::endl;
    std::cout << "    k_max: " << std::scientific << k_max << std::defaultfloat << std::endl;
    std::cout << "    k_mean: " << std::scientific << k_mean << std::defaultfloat << std::endl;
    std::cout << "    ε_max: " << std::scientific << eps_max << std::defaultfloat << std::endl;
    std::cout << "    νt_max: " << std::scientific << nu_t_max << std::defaultfloat 
              << " (ratio to ν: " << nu_t_max/nu << "x)" << std::endl;
    std::cout << "    νt_mean: " << std::scientific << nu_t_mean << std::defaultfloat << std::endl;
    
    // Wall y+ diagnostics
    auto [max_yplus, mean_yplus] = solver::ResidualMonitor::wallYPlus(u, v, mesh, nu);
    std::cout << "  Wall y+ statistics:" << std::endl;
    std::cout << "    y+_max: " << std::fixed << std::setprecision(2) << max_yplus << std::endl;
    std::cout << "    y+_mean: " << std::fixed << std::setprecision(2) << mean_yplus << std::endl;
    if (max_yplus > 11.25) {
        std::cout << "    (in log-layer region - wall functions appropriate)" << std::endl;
    } else {
        std::cout << "    (in viscous sublayer - low-Re treatment required)" << std::endl;
    }

    // Save fields to CSV for analysis
    io::CSVWriter::writeFieldToCSV(csvOutputDir + "/velocity_u.csv", u_cc.raw());
    io::CSVWriter::writeFieldToCSV(csvOutputDir + "/velocity_v.csv", v_cc.raw());
    io::CSVWriter::writeFieldToCSV(csvOutputDir + "/pressure.csv", p.raw());
    io::CSVWriter::writeFieldToCSV(csvOutputDir + "/turbulent_k.csv", k_field.raw());
    io::CSVWriter::writeFieldToCSV(csvOutputDir + "/turbulent_epsilon.csv", eps_field.raw());
    io::CSVWriter::writeFieldToCSV(csvOutputDir + "/turbulent_viscosity.csv", nu_t_field.raw());

    // Write centerline profiles
    ::io::BenchmarkOutput::writeVerticalCenterlineU(csvOutputDir + "/centerline_u.csv", mesh, u_cc);
    ::io::BenchmarkOutput::writeHorizontalCenterlineV(csvOutputDir + "/centerline_v.csv", mesh, v_cc);

    // Finalize animation
    if (cfg.useAnimation) {
        pvdWriter.finalize();
        pvdWriter.writePVD();
    }

    std::cout << "  Output files: " << csvOutputDir << "/" << std::endl;
    std::cout << "=== High Reynolds Cavity Complete ===" << std::endl << std::endl;
}

}  // namespace cases
