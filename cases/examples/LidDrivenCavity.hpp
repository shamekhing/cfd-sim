#pragma once

#include <iostream>
#include <cmath>
#include <filesystem>
#include <chrono>
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
#include "../../include/solver/VortexMetrics.hpp"
#include "../../include/physics/TurbulenceKEpsilon.hpp"
#include "../../include/io/VTKWriter.hpp"
#include "../../include/io/ColorMap.hpp"
#include "../../include/io/PVDWriter.hpp"
#include "../../include/io/VTSReader.hpp"
#include "../../include/io/BenchmarkOutput.hpp"
#include "../../include/numerics/VortexAnalytics.hpp"
#include "../../tests/validation/GhiaReference.hpp"
#include "../../tests/validation/LidCavityVortexRef.hpp"

namespace cases {

inline void runLidDrivenCavity(
    const config::LidDrivenCavityConfig& cfg,
    const std::string& outputDir = "build/output",
    const std::string& name = "output")
{
    const int Nx = cfg.Nx;
    const int Ny = cfg.Ny;
    const double Re = cfg.Re;
    const double nu = cfg.nu();
    const std::string outputRoot = outputDir;
    const std::string runOutputDir = outputRoot + "/" + name;
    const std::string csvOutputDir = runOutputDir + "/csv";
    const std::string vtsOutputDir = runOutputDir + "/vts";
    const bool useAnimation = cfg.useAnimation;

    core::Mesh mesh(Nx, Ny, cfg.Lx, cfg.Ly);

    core::UField u(mesh);
    core::VField v(mesh);
    core::Field p(mesh);

    u.setZero();
    v.setZero();
    p.setZero();
    solver::StaggeredSIMPLE::applyVelocityBC(u, v, Nx, Ny);

    // Initialize turbulence model if enabled
    physics::TurbulenceKEpsilon turbulence;
    physics::TurbulenceKEpsilon* turb_ptr = nullptr;
    if (cfg.enable_turbulence) {
        turbulence.initialize(mesh, config::TurbulenceConfig::standardKEpsilon());
        turb_ptr = &turbulence;
        if (cfg.verbosity >= 1) {
            std::cout << "Turbulence k-epsilon model enabled\n";
        }
    }

    if (!useAnimation) {
        // Steady-state mode - run to convergence
        solver::StaggeredSIMPLE::solveSteady(mesh, u, v, p, nu, cfg.maxSteadyIter,
            cfg.alpha_u, cfg.alpha_v, cfg.alpha_p, cfg.verbosity, turb_ptr);
    } else {
        // Transient animation mode - time-stepping
        const double dt = cfg.dt;
        const int numTimeSteps = cfg.numTimeSteps;
        const int itersPerTimeStep = cfg.itersPerTimeStep;

        // Initialize output directories and animation writer
        std::filesystem::create_directories(runOutputDir);
        std::filesystem::create_directories(csvOutputDir);
        std::filesystem::create_directories(vtsOutputDir);
        io::PVDWriter pvdWriter(vtsOutputDir, name, cfg.verbosity, runOutputDir, cfg.colorScheme, cfg.useColorMap);
        
        // Check if resuming from existing frames
        int startFrame = pvdWriter.getStartFrame();
        std::string lastFrameFile = pvdWriter.getLastFrameFilename();
        
        // Load field state from last frame if resuming
        if (startFrame > 0 && !lastFrameFile.empty()) {
            if (cfg.verbosity > 0) {
                std::cout << "  Loading field state from: " << lastFrameFile << std::endl;
            }
            
            // Load cell-centered data from VTS file
            core::Field u_cc(mesh);
            core::Field v_cc(mesh);
            core::Field p_loaded(mesh);
            
            if (io::VTSReader::read(lastFrameFile, mesh, u_cc, v_cc, p_loaded)) {
                // Copy loaded cell-centered values to staggered fields
                // (Simplified: use cell-centered values directly - proper interpolation could be added)
                for (int j = 0; j < Ny; ++j) {
                    for (int i = 0; i < Nx; ++i) {
                        u(i, j) = u_cc(i, j);
                        v(i, j) = v_cc(i, j);
                        p(i, j) = p_loaded(i, j);
                    }
                }
                // Reapply boundary conditions
                solver::StaggeredSIMPLE::applyVelocityBC(u, v, Nx, Ny);
                
                if (cfg.verbosity > 0) {
                    std::cout << "  Successfully loaded field state from frame " << (startFrame - 1) << std::endl;
                }
            } else {
                std::cerr << "  Warning: Failed to load field state, starting from zero" << std::endl;
                startFrame = 0;
            }
        }

        core::UField u_old(mesh);
        core::VField v_old(mesh);

        auto anim_start = std::chrono::high_resolution_clock::now();

        for (int tstep = startFrame; tstep < numTimeSteps; ++tstep) {
            u_old.raw() = u.raw();
            v_old.raw() = v.raw();
            
            // Solve transient silently (verbosity=0) - animation progress shown separately
            solver::StaggeredSIMPLE::solveTransient(mesh, u, v, p, u_old, v_old, nu, dt,
                itersPerTimeStep, cfg.alpha_u, cfg.alpha_v, cfg.alpha_p, 0, turb_ptr);
            
            // Output frame at each timestep
            core::Field u_cc = u.cellCenterView(mesh);
            core::Field v_cc = v.cellCenterView(mesh);
            double time = (tstep + 1) * dt;
            pvdWriter.writeTimestep(tstep, time, mesh, u_cc, v_cc, p);

            // Print single-line progress
            auto current_time = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(current_time - anim_start).count();
            double avg_time = elapsed / (tstep - startFrame + 1);
            int remaining = numTimeSteps - tstep - 1;
            double eta = avg_time * remaining;
            
            int progress = ((tstep + 1) * 100) / numTimeSteps;
            std::string prefix = (startFrame > 0) ? "  Resuming: " : "  Animation: ";
            std::cout << "\r" << prefix << (tstep + 1) << "/" << numTimeSteps 
                      << " [" << std::string(progress / 5, '#') 
                      << std::string(20 - progress / 5, '-') << "] " 
                      << progress << "% | "
                      << std::fixed << std::setprecision(1)
                      << "T:" << elapsed << "s | ETA:" << eta << "s"
                      << std::flush;
        }
        std::cout << std::endl;

        // Wait for all async writes to complete
        pvdWriter.finalize();
        
        // Write the PVD collection file
        pvdWriter.writePVD();
    }

    double massRes = solver::ResidualMonitor::massResidualL2Staggered(mesh, u, v);

    if (cfg.verbosity >= 1) {
        std::cout << "Lid Driven Cavity Completed (staggered)\n";
        std::cout << "Re = " << Re << "\n";
        std::cout << "Mass residual L2 (staggered) = " << massRes << "\n";
    }

    core::Field u_cc = u.cellCenterView(mesh);
    core::Field v_cc = v.cellCenterView(mesh);

    double maxVelInterior = 0.0;
    for (int j = 1; j < Ny - 1; ++j) {
        for (int i = 1; i < Nx - 1; ++i) {
            double mag = std::sqrt(u_cc(i, j) * u_cc(i, j) + v_cc(i, j) * v_cc(i, j));
            if (mag > maxVelInterior) maxVelInterior = mag;
        }
    }
    if (cfg.verbosity >= 1) {
        std::cout << "Max velocity magnitude (interior, cell-centered) = " << maxVelInterior << std::endl;
    }

    if (cfg.verbosity >= 1) {
        bool has_reference = false;
        double errU = 0.0;
        double errV = 0.0;
        if (std::abs(Re - 100.0) < 1e-6) {
            errU = validation::centerlineErrorU(mesh, u_cc, validation::GhiaRe100::y_ref, validation::GhiaRe100::u_ref);
            errV = validation::centerlineErrorV(mesh, v_cc, validation::GhiaRe100::x_ref, validation::GhiaRe100::v_ref);
            has_reference = true;
            std::cout << "Ghia L2 error (Re=100, u centerline) = " << errU << std::endl;
            std::cout << "Ghia L2 error (Re=100, v centerline) = " << errV << std::endl;
        } else if (std::abs(Re - 1000.0) < 1e-6) {
            errU = validation::centerlineErrorU(mesh, u_cc, validation::GhiaRe1000::y_ref, validation::GhiaRe1000::u_ref);
            errV = validation::centerlineErrorV(mesh, v_cc, validation::GhiaRe1000::x_ref, validation::GhiaRe1000::v_ref);
            has_reference = true;
            std::cout << "Ghia L2 error (Re=1000, u centerline) = " << errU << std::endl;
            std::cout << "Ghia L2 error (Re=1000, v centerline) = " << errV << std::endl;
        }

        if (cfg.validationMode && has_reference) {
            const bool mass_ok = massRes < 1e-4;
            const bool u_ok = errU < 0.20;
            const bool v_ok = errV < 0.20;
            std::cout << "Validation status: "
                      << ((mass_ok && u_ok && v_ok) ? "PASS" : "FAIL")
                      << " (mass<1e-4, uL2<0.20, vL2<0.20)" << std::endl;
        }
        std::cout << "VTK file written: " << runOutputDir << "/" << name << ".vtk\n";
    }

    std::filesystem::create_directories(runOutputDir);
    std::filesystem::create_directories(csvOutputDir);

    // Standard VTK output
    if (cfg.useColorMap) {
        io::VTKWriter::writeWithColorMap(runOutputDir + "/" + name + ".vtk", mesh, u_cc, v_cc, p, cfg.colorScheme);
    } else {
        io::VTKWriter::write(runOutputDir + "/" + name + ".vtk", mesh, u_cc, v_cc, p);
    }
    
    io::BenchmarkOutput::writeVerticalCenterlineU(csvOutputDir + "/centerline_u.csv", mesh, u_cc);
    io::BenchmarkOutput::writeHorizontalCenterlineV(csvOutputDir + "/centerline_v.csv", mesh, v_cc);
    
    // Vortex metrics calculation
    if (cfg.computeVortexMetrics) {
        if (cfg.verbosity >= 1) {
            std::cout << "\nComputing vortex metrics...\n";
        }
        
        // Compute vorticity and streamfunction
        core::Field omega = numerics::VortexAnalytics::computeVorticity(mesh, u_cc, v_cc);
        core::Field psi = numerics::VortexAnalytics::computeStreamfunction(mesh, u_cc, v_cc);
        
        // Find primary vortex center
        auto [vortex_x, vortex_y] = numerics::VortexAnalytics::findPrimaryVortexCenter(mesh, p);
        
        // Compute circulation and enstrophy
        double circulation = solver::VortexMetrics::computeCirculation(mesh, omega, vortex_x, vortex_y, 0.15);
        double enstrophy = solver::VortexMetrics::computeEnstrophy(mesh, omega);
        double max_vort = solver::VortexMetrics::computeMaxVorticity(omega);
        
        // Find secondary vortices
        auto secondary_vortices = numerics::VortexAnalytics::findSecondaryVortices(
            mesh, omega, cfg.vortexDetectionThreshold);
        
        if (cfg.verbosity >= 1) {
            std::cout << "  Primary vortex center: (" << vortex_x << ", " << vortex_y << ")\n";
            std::cout << "  Circulation: " << circulation << "\n";
            std::cout << "  Enstrophy: " << enstrophy << "\n";
            std::cout << "  Max vorticity: " << max_vort << "\n";
            std::cout << "  Secondary vortices found: " << secondary_vortices.size() << "\n";
            for (size_t i = 0; i < secondary_vortices.size(); ++i) {
                std::cout << "    Vortex " << (i+1) << ": (" 
                          << secondary_vortices[i].first << ", "
                          << secondary_vortices[i].second << ")\n";
            }
        }
        
        // Validate against reference
        auto [ref_x, ref_y] = validation::getVortexReferenceCenter(Re);
        auto vortex_result = solver::VortexMetrics::validateVortexMetrics(
            vortex_x, vortex_y, ref_x, ref_y, circulation, 0.0);
        
        if (cfg.verbosity >= 1) {
            std::cout << "  Vortex center error vs reference: " << vortex_result.center_error << "\n";
        }
        
        // Output vortex fields
        if (cfg.outputVorticity || cfg.outputStreamfunction) {
            io::VTKWriter::writeVortexFields(
                runOutputDir + "/" + name + "_vortex.vtk",
                mesh, omega, psi, cfg.colorScheme);
            
            if (cfg.verbosity >= 1) {
                std::cout << "  Vortex fields written to: " << name << "_vortex.vtk\n";
            }
        }
        
        // Write vortex metrics to CSV
        std::ofstream vortex_csv(csvOutputDir + "/vortex_metrics.csv");
        vortex_csv << "metric,value,reference\n";
        vortex_csv << "primary_vortex_x," << vortex_x << "," << ref_x << "\n";
        vortex_csv << "primary_vortex_y," << vortex_y << "," << ref_y << "\n";
        vortex_csv << "circulation," << circulation << ",0.0\n";
        vortex_csv << "enstrophy," << enstrophy << ",0.0\n";
        vortex_csv << "max_vorticity," << max_vort << ",0.0\n";
        vortex_csv << "center_error," << vortex_result.center_error << ",0.0\n";
        vortex_csv.close();
    }
    
    // Output turbulence fields if enabled
    if (cfg.enable_turbulence && turb_ptr) {
        if (cfg.verbosity >= 1) {
            std::cout << "Writing turbulence fields to VTK...\n";
        }
        const auto& k_field = turb_ptr->k();
        const auto& eps_field = turb_ptr->epsilon();
        const auto& nu_t_field = turb_ptr->turbulentViscosity();
        
        io::VTKWriter::writeWithTurbulence(
            runOutputDir + "/" + name + "_turbulence.vtk",
            mesh,
            u_cc, v_cc, p,
            k_field, eps_field, nu_t_field);
        
        if (cfg.verbosity >= 1) {
            // Print turbulence statistics
            double k_max = 0.0, k_mean = 0.0;
            double eps_max = 0.0;
            for (int idx = 0; idx < k_field.size(); ++idx) {
                k_max = std::max(k_max, k_field[idx]);
                k_mean += k_field[idx];
                eps_max = std::max(eps_max, eps_field[idx]);
            }
            k_mean /= k_field.size();
            std::cout << "Turbulence statistics:\n";
            std::cout << "  k_max: " << std::scientific << k_max << std::defaultfloat << "\n";
            std::cout << "  k_mean: " << std::scientific << k_mean << std::defaultfloat << "\n";
            std::cout << "  ε_max: " << std::scientific << eps_max << std::defaultfloat << "\n";
        }
    }
}

} // namespace cases
