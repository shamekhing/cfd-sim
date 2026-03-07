#include <iostream>
#include <cmath>
#include <cassert>

#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"
#include "../../include/core/Indexing.hpp"
#include "../../include/assembly/PressureAssembler.hpp"
#include "../../include/linear/BiCGSTAB.hpp"

/**
 * Unit test for Poisson solver using Method of Manufactured Solutions (MMS)
 * 
 * Tests the Laplacian operator and pressure assembler.
 * 
 * Manufactured solution:
 *   φ(x,y) = sin(πx) sin(πy)
 *   ∇²φ = -2π² sin(πx) sin(πy)
 * 
 * Boundary conditions: φ = 0 on all boundaries (homogeneous Dirichlet)
 * 
 * Acceptance criteria: L2 error < 0.01
 */

int main() {
    std::cout << "=============================================================\n";
    std::cout << "  Unit Test: Poisson Solver (MMS Verification)\n";
    std::cout << "=============================================================\n\n";

    const int Nx = 64;
    const int Ny = 64;

    core::Mesh mesh(Nx, Ny, 1.0, 1.0);
    core::Field rhs(mesh);

    std::cout << "Configuration:\n";
    std::cout << "  Grid: " << Nx << " x " << Ny << "\n";
    std::cout << "  Domain: [0,1] x [0,1]\n";
    std::cout << "  Manufactured solution: φ = sin(πx) sin(πy)\n";
    std::cout << "  Laplacian: ∇²φ = -2π² sin(πx) sin(πy)\n\n";

    // Assemble RHS from Laplacian of manufactured solution
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double x = mesh.x(i);
            double y = mesh.y(j);

            rhs(i, j) = -2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
        }
    }

    std::cout << "Assembling linear system...\n";
    auto system = assembly::PressureAssembler::assemble(mesh, rhs, 0.0);

    std::cout << "Solving with BiCGSTAB...\n";
    linear::Vector phi = linear::BiCGSTABSolver().solve(system.A, system.b);

    std::cout << "Computing error...\n";
    double error_l2 = 0.0;
    double error_max = 0.0;
    int count = 0;

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double x = mesh.x(i);
            double y = mesh.y(j);
            double exact = std::sin(M_PI * x) * std::sin(M_PI * y);
            
            int P = core::id(i, j, Nx);
            double error = phi(P) - exact;
            
            error_l2 += error * error;
            error_max = std::max(error_max, std::abs(error));
            ++count;
        }
    }

    error_l2 = std::sqrt(error_l2 / count);

    std::cout << "\nResults:\n";
    std::cout << "  L2 error:  " << error_l2 << "\n";
    std::cout << "  Max error: " << error_max << "\n\n";

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
