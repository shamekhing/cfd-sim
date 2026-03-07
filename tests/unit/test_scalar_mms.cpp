#include <iostream>
#include <cmath>
#include <iomanip>

#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"
#include "../../include/assembly/GenericScalarAssembler.hpp"
#include "../../include/linear/BiCGSTAB.hpp"

/**
 * Unit test for scalar transport using Method of Manufactured Solutions (MMS)
 * 
 * Tests the scalar transport equation:
 *   ∂T/∂t = ∇·(α∇T)
 * 
 * Manufactured solution:
 *   T(x,y,t) = cos(2πx) cos(2πy) exp(-4π²αt)
 * 
 * Acceptance criteria: L2 error < 0.01
 */

int main() {
    std::cout << "=============================================================\n";
    std::cout << "  Unit Test: Scalar Transport (MMS Verification)\n";
    std::cout << "=============================================================\n\n";

    const int Nx = 64;
    const int Ny = 64;
    const double Lx = 1.0;
    const double Ly = 1.0;

    core::Mesh mesh(Nx, Ny, Lx, Ly);
    
    core::Field T(mesh);
    core::Field T_old(mesh);
    
    double alpha = 0.01;  // Thermal diffusivity
    
    // Zero velocity field
    core::Field u(mesh);
    core::Field v(mesh);
    u.setZero();
    v.setZero();
    
    // Initialize with manufactured solution at t=0
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double x = mesh.x(i);
            double y = mesh.y(j);
            T(i, j) = std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
        }
    }
    T_old = T;
    
    double dt = 1e-4;
    int num_steps = 100;
    double t = 0.0;
    double t_final = num_steps * dt;
    double decay_factor = std::exp(-4.0 * M_PI * M_PI * alpha * t_final);
    
    std::cout << "Configuration:\n";
    std::cout << "  Grid: " << Nx << " x " << Ny << "\n";
    std::cout << "  Domain: [0,1] x [0,1]\n";
    std::cout << "  Thermal diffusivity α: " << alpha << "\n";
    std::cout << "  Time step: " << dt << "\n";
    std::cout << "  Number of steps: " << num_steps << "\n";
    std::cout << "  Total time: " << t_final << "\n";
    std::cout << "  Expected decay factor: " << decay_factor << "\n\n";
    
    std::cout << "Running time integration...\n";
    
    for (int step = 0; step < num_steps; ++step) {
        T_old = T;
        
        auto system = assembly::GenericScalarAssembler::assemble(
            mesh, T, u, v, alpha, dt, T_old,
            assembly::ConvectionScheme::Upwind,
            nullptr, true, 0.01);
        
        linear::Vector T_new = linear::BiCGSTABSolver().solve(system.A, system.b);
        T.raw() = T_new;
        
        t += dt;
        
        if ((step + 1) % 20 == 0) {
            std::cout << "  Step " << (step + 1) << " / " << num_steps
                      << ", t = " << std::setprecision(6) << t << "\n";
        }
    }
    
    std::cout << "\nComputing error...\n";
    
    double error_l2 = 0.0;
    double error_max = 0.0;
    double T_max_exact = 0.0;
    int count = 0;
    
    for (int j = 1; j < Ny - 1; ++j) {
        for (int i = 1; i < Nx - 1; ++i) {
            double x = mesh.x(i);
            double y = mesh.y(j);
            
            double T_exact = decay_factor * std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
            double T_num = T(i, j);
            double error = T_num - T_exact;
            
            error_l2 += error * error;
            error_max = std::max(error_max, std::abs(error));
            T_max_exact = std::max(T_max_exact, std::abs(T_exact));
            ++count;
        }
    }
    
    error_l2 = std::sqrt(error_l2 / count);
    double relative_error = (error_l2 / T_max_exact) * 100.0;
    
    std::cout << "\nResults:\n";
    std::cout << "  L2 error: " << std::setprecision(8) << error_l2 << "\n";
    std::cout << "  Max error: " << error_max << "\n";
    std::cout << "  Relative error: " << relative_error << " %\n\n";
    
    // Acceptance criteria
    const double threshold = 0.01;
    bool pass = error_l2 < threshold;
    
    std::cout << "Acceptance Criteria:\n";
    std::cout << "  [" << (pass ? "✓" : "✗") << "] L2 error < " << threshold 
              << " (" << error_l2 << ")\n\n";
    
    std::cout << "=============================================================\n";
    std::cout << "  Test Result: " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "=============================================================\n";
    
    return pass ? 0 : 1;
}
