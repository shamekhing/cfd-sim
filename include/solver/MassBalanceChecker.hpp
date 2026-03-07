#pragma once

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "../include/core/UField.hpp"
#include "../include/core/VField.hpp"
#include "../include/linear/Vector.hpp"

namespace solver {

/**
 * Mass balance and continuity checking utilities
 * 
 * Verifies that the incompressibility constraint (∇·u = 0) is satisfied
 * at the discrete level. Mass balance violations indicate:
 *  - Convergence issues
 *  - Boundary condition problems
 *  - Numerical errors in pressure-velocity coupling
 * 
 * For staggered grids, checks:
 *  - Cell-wise mass balance: (u_w - u_e)*dy + (v_s - v_n)*dx ≈ 0
 *  - Global total: Σ mass = should be near zero
 *  - Violation count: how many cells violate continuity
 * 
 * For collocated grids, checks:
 *  - Divergence: div(u,v) = ∂u/∂x + ∂v/∂y ≈ 0
 */
class MassBalanceChecker {
public:
    /**
     * Result of a mass balance check
     */
    struct MassBalanceResult {
        double globalMassError;      // Σ |mass_i| over entire domain
        double maxCellMassError;     // max |mass_i| in any single cell
        double avgCellMassError;     // average |mass_i| over domain
        int violationCount;          // number of cells violating continuity
        double violationFraction;    // violationCount / totalCells
        int totalCells;              // total interior cells checked
    };

    /**
     * Check staggered grid mass balance
     * 
     * For each interior cell (i,j), compute continuity error:
     *   mass_error = (u_west - u_east) * dy + (v_south - v_north) * dx
     * 
     * This should be near zero for incompressible flow.
     * 
     * @param mesh   Structured grid
     * @param u      U-velocity at face centers
     * @param v      V-velocity at face centers
     * @param threshold  Cells with |mass| > threshold count as violations
     * @return  Statistics on mass balance
     */
    static MassBalanceResult checkStaggeredMassBalance(
        const core::Mesh& mesh,
        const core::UField& u,
        const core::VField& v,
        linear::Scalar threshold = 1e-6)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        MassBalanceResult result;
        result.globalMassError = 0.0;
        result.maxCellMassError = 0.0;
        result.avgCellMassError = 0.0;
        result.violationCount = 0;
        result.totalCells = (Nx - 2) * (Ny - 2);

        if (result.totalCells == 0) {
            result.violationFraction = 0.0;
            return result;
        }

        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                // Get velocities at cell faces
                // West face: u(i, j)
                // East face: u(i+1, j)
                // South face: v(i, j)
                // North face: v(i, j+1)
                
                linear::Scalar u_w = u(i, j);
                linear::Scalar u_e = u(i + 1, j);
                linear::Scalar v_s = v(i, j);
                linear::Scalar v_n = v(i, j + 1);

                // Mass conservation for this cell (should be ~0):
                // ∇·u = ∂u/∂x + ∂v/∂y ≈ (u_e - u_w)/dx + (v_n - v_s)/dy
                // Integrated over cell: (u_e - u_w)*dy + (v_n - v_s)*dx ≈ 0
                linear::Scalar cellMass = (u_e - u_w) * dy + (v_n - v_s) * dx;
                linear::Scalar celMassAbs = std::abs(cellMass);

                result.globalMassError += celMassAbs;
                result.maxCellMassError = std::max(result.maxCellMassError, celMassAbs);
                result.avgCellMassError += celMassAbs;

                if (celMassAbs > threshold) {
                    result.violationCount++;
                }
            }
        }

        result.avgCellMassError /= result.totalCells;
        result.violationFraction = static_cast<double>(result.violationCount) / result.totalCells;

        return result;
    }

    /**
     * Check collocated grid divergence (for debugging/reference)
     * 
     * Computes divergence using central differences:
     *   div = (u_e - u_w)/(2*dx) + (v_n - v_s)/(2*dy)
     * 
     * Returns L2 norm of divergence field.
     * 
     * @param mesh   Structured grid
     * @param u, v   Velocity components (cell-centered for this version)
     * @return  L2 norm of divergence
     */
    static linear::Scalar checkCollocatedDivergence(
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        linear::Scalar divergenceL2 = 0.0;
        int count = 0;

        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                linear::Scalar du_dx = (u(i + 1, j) - u(i - 1, j)) / (2.0 * dx);
                linear::Scalar dv_dy = (v(i, j + 1) - v(i, j - 1)) / (2.0 * dy);
                linear::Scalar div = du_dx + dv_dy;

                divergenceL2 += div * div;
                ++count;
            }
        }

        return (count > 0) ? std::sqrt(divergenceL2 / count) : 0.0;
    }

    /**
     * Print formatted mass balance report
     * 
     * @param result    Result from checkStaggeredMassBalance
     * @param threshold Used in check (for reporting)
     */
    static void printReport(
        const MassBalanceResult& result,
        linear::Scalar threshold = 1e-6)
    {
        std::cout << "\n  ===== Mass Balance Report =====" << std::endl;
        std::cout << "  Total interior cells: " << result.totalCells << std::endl;
        std::cout << "  Global mass error (L1): " << std::scientific << std::setprecision(3)
                  << result.globalMassError << std::defaultfloat << std::endl;
        std::cout << "  Max cell mass error: " << std::scientific << std::setprecision(3)
                  << result.maxCellMassError << std::defaultfloat << std::endl;
        std::cout << "  Avg cell mass error: " << std::scientific << std::setprecision(3)
                  << result.avgCellMassError << std::defaultfloat << std::endl;
        std::cout << "  Continuity violations (|mass| > " << threshold << "): " 
                  << result.violationCount << " / " << result.totalCells
                  << " (" << std::fixed << std::setprecision(1) 
                  << (result.violationFraction * 100.0) << "%)" << std::endl;
        
        // Status indicator
        std::string status = "GOOD";
        if (result.violationFraction > 0.5) status = "POOR";
        else if (result.violationFraction > 0.1) status = "FAIR";
        
        std::cout << "  Status: " << status << std::endl;
        std::cout << "  ==============================\n" << std::endl;
    }

    /**
     * Detailed per-cell violation report (warning: verbose!)
     * 
     * Lists all cells violating threshold, with location and error magnitude
     * 
     * @param mesh       Structured grid
     * @param u, v       Velocity fields
     * @param threshold  Report cells with |error| > threshold
     * @param maxLines   Max violations to print (default unlimited)
     */
    static void printViolations(
        const core::Mesh& mesh,
        const core::UField& u,
        const core::VField& v,
        linear::Scalar threshold = 1e-6,
        int maxLines = 50)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        int count = 0;
        std::cout << "\n  Mass Balance Violations (|error| > " << threshold << "):\n";
        std::cout << "  Cell (i,j)  |  Mass Error\n";
        std::cout << "  " << std::string(30, '-') << std::endl;

        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                linear::Scalar u_w = u(i, j);
                linear::Scalar u_e = u(i + 1, j);
                linear::Scalar v_s = v(i, j);
                linear::Scalar v_n = v(i, j + 1);

                linear::Scalar cellMass = std::abs((u_e - u_w) * dy + (v_n - v_s) * dx);

                if (cellMass > threshold) {
                    std::cout << "  (" << std::setw(3) << i << "," << std::setw(3) << j << ")  |  "
                              << std::scientific << cellMass << std::defaultfloat << std::endl;
                    count++;
                    if (count >= maxLines) {
                        std::cout << "  ... (showing first " << maxLines << " violations) ...\n" << std::endl;
                        return;
                    }
                }
            }
        }
        std::cout << std::endl;
    }

};

}  // namespace solver