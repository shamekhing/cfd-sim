#include <iostream>
#include <cassert>
#include <cmath>
#include <filesystem>

#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"
#include "../../include/core/UField.hpp"
#include "../../include/core/VField.hpp"
#include "../../include/config/CaseConfigs.hpp"
#include "../../include/solver/StaggeredSIMPLE.hpp"
#include "../../include/solver/ResidualMonitor.hpp"
#include "../../include/io/VTKWriter.hpp"
#include "../../include/io/BenchmarkOutput.hpp"
#include "GhiaReference.hpp"

/**
 * Automated validation test for Ghia Re=100 lid-driven cavity benchmark.
 * 
 * Reference: Ghia et al. (1982) JCP 48:387-411
 * DOI: 10.1016/0021-9991(82)90058-4
 * 
 * Acceptance criteria:
 *   - L2 error on vertical centerline u(y): < 0.02
 *   - L2 error on horizontal centerline v(x): < 0.02
 *   - Mass conservation residual: < 1e-4
 */

int main() {
    std::cout << "=============================================================\n";
    std::cout << "  Validation Test: Lid-Driven Cavity Re=100 (Ghia 1982)\n";
    std::cout << "=============================================================\n\n";

    // Configuration for Ghia Re=100 benchmark
    // Grid: 129x129 (minimum for convergence per Ghia paper)
    config::LidDrivenCavityConfig cfg;
    cfg.Re = 100.0;
    cfg.Nx = 129;
    cfg.Ny = 129;
    cfg.Lx = 1.0;
    cfg.Ly = 1.0;
    cfg.maxSteadyIter = 10000;
    cfg.alpha_u = 0.7;
    cfg.alpha_v = 0.7;
    cfg.alpha_p = 0.3;
    cfg.tolerance = 1e-6;
    cfg.useAnimation = false;
    cfg.validationMode = true;
    cfg.enable_turbulence = false;
    cfg.verbosity = 1;

    const double nu = cfg.nu();
    
    std::cout << "Configuration:\n";
    std::cout << "  Grid: " << cfg.Nx << " x " << cfg.Ny << "\n";
    std::cout << "  Reynolds number: " << cfg.Re << "\n";
    std::cout << "  Viscosity: " << nu << "\n";
    std::cout << "  Max iterations: " << cfg.maxSteadyIter << "\n";
    std::cout << "  Convergence tolerance: " << cfg.tolerance << "\n\n";

    // Setup mesh and fields
    core::Mesh mesh(cfg.Nx, cfg.Ny, cfg.Lx, cfg.Ly);
    core::UField u(mesh);
    core::VField v(mesh);
    core::Field p(mesh);

    u.setZero();
    v.setZero();
    p.setZero();

    std::cout << "Running SIMPLE solver...\n\n";
    
    // Use the solveSteady method for automated convergence
    solver::StaggeredSIMPLE::solveSteady(
        mesh, u, v, p, nu, cfg.maxSteadyIter,
        cfg.alpha_u, cfg.alpha_v, cfg.alpha_p,
        cfg.verbosity, nullptr);

    // Convert to cell-centered for comparison
    core::Field u_cc(mesh), v_cc(mesh);
    for (int j = 0; j < cfg.Ny; ++j) {
        for (int i = 0; i < cfg.Nx; ++i) {
            u_cc(i, j) = 0.5 * (u(i, j) + u(i + 1, j));
            v_cc(i, j) = 0.5 * (v(i, j) + v(i, j + 1));
        }
    }

    // Compare with Ghia benchmark data
    std::cout << "Comparing with Ghia Re=100 reference data...\n";
    
    const double errU = validation::centerlineErrorU(
        mesh, u_cc, validation::GhiaRe100::y_ref, validation::GhiaRe100::u_ref);
    
    const double errV = validation::centerlineErrorV(
        mesh, v_cc, validation::GhiaRe100::x_ref, validation::GhiaRe100::v_ref);

    const double massRes = solver::ResidualMonitor::massResidualL2Staggered(mesh, u, v);

    std::cout << "\nValidation Metrics:\n";
    std::cout << "  L2 error (u centerline): " << errU << "\n";
    std::cout << "  L2 error (v centerline): " << errV << "\n";
    std::cout << "  Mass residual:           " << massRes << "\n\n";

    // Acceptance criteria thresholds
    const double u_threshold = 0.02;  // 2% L2 error
    const double v_threshold = 0.02;  // 2% L2 error
    const double mass_threshold = 1e-4;

    bool u_pass = errU < u_threshold;
    bool v_pass = errV < v_threshold;
    bool mass_pass = massRes < mass_threshold;
    bool all_pass = u_pass && v_pass && mass_pass;

    std::cout << "Acceptance Criteria:\n";
    std::cout << "  [" << (u_pass ? "✓" : "✗") << "] u centerline L2 < " << u_threshold 
              << " (" << errU << ")\n";
    std::cout << "  [" << (v_pass ? "✓" : "✗") << "] v centerline L2 < " << v_threshold 
              << " (" << errV << ")\n";
    std::cout << "  [" << (mass_pass ? "✓" : "✗") << "] mass residual < " << mass_threshold 
              << " (" << massRes << ")\n\n";

    std::cout << "=============================================================\n";
    std::cout << "  Test Result: " << (all_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "=============================================================\n";

    // Write output for inspection
    std::filesystem::create_directories("build/output/validation");
    io::VTKWriter::write("build/output/validation/cavity_re100.vtk", mesh, u_cc, v_cc, p);
    io::BenchmarkOutput::writeVerticalCenterlineU("build/output/validation/cavity_re100_u.csv", mesh, u_cc);
    io::BenchmarkOutput::writeHorizontalCenterlineV("build/output/validation/cavity_re100_v.csv", mesh, v_cc);
    
    std::cout << "\nOutput written to build/output/validation/cavity_re100.*\n";

    // Return non-zero exit code if test fails
    return all_pass ? 0 : 1;
}
