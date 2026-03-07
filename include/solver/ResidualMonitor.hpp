#pragma once

#include <cmath>
#include <stdexcept>
#include <utility>
#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "../include/core/UField.hpp"
#include "../include/core/VField.hpp"
#include "../include/core/Indexing.hpp"
#include "../include/linear/Vector.hpp"

namespace solver {

class ResidualMonitor {
public:
    // Staggered: mass imbalance per cell (u_w - u_e)*dy + (v_s - v_n)*dx; L2 over interior
    static linear::Scalar massResidualL2Staggered(const core::Mesh& mesh,
                                                const core::UField& u,
                                                const core::VField& v)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();
        linear::Scalar sum = 0.0;
        int count = 0;
        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                linear::Scalar u_w = u(i, j);
                linear::Scalar u_e = u(i + 1, j);
                linear::Scalar v_s = v(i, j);
                linear::Scalar v_n = v(i, j + 1);
                linear::Scalar mass = (u_w - u_e) * dy + (v_s - v_n) * dx;
                sum += mass * mass;
                ++count;
            }
        }
        return (count > 0) ? std::sqrt(sum / count) : 0.0;
    }
    // Compute residual norm ||Ax - b|| for a linear system
    // Uses SparseMatrix::operator* to multiply A*x, then computes ||A*x - b||_2
    static linear::Scalar residualNorm(const linear::SparseMatrix& A,
                                     const linear::Vector& x,
                                     const linear::Vector& b)
    {
        if (x.size() != b.size() || x.size() != A.rows()) {
            return 0.0;
        }
        
        // Compute Ax using sparse matrix multiplication
        linear::Vector Ax = A * x;
        
        // Compute residual: r = Ax - b
        linear::Vector residual = Ax;
        residual.raw() -= b.raw();
        
        // Return L2 norm: ||r||_2 = sqrt(r · r)
        linear::Scalar norm_sq = 0.0;
        for (int i = 0; i < residual.size(); ++i) {
            norm_sq += residual[i] * residual[i];
        }
        
        return std::sqrt(norm_sq);
    }

    // Compute divergence L2 norm over interior cells
    static linear::Scalar divergenceL2(const core::Mesh& mesh,
                                     const core::Field& u,
                                     const core::Field& v)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        linear::Scalar sum = 0.0;
        int count = 0;

        for (int j = 1; j < Ny-1; ++j)
        {
            for (int i = 1; i < Nx-1; ++i)
            {
                                linear::Scalar div =
                    (u(i+1,j) - u(i-1,j)) / (2.0*dx)
                  + (v(i,j+1) - v(i,j-1)) / (2.0*dy);
                sum += div * div;
                ++count;
            }
        }

        return (count > 0) ? std::sqrt(sum / count) : 0.0;
    }

    // Pressure correction residual: ||A p' - b|| (for convergence logging)
    static linear::Scalar pressureResidualNorm(const linear::SparseMatrix& A,
                                             const linear::Vector& pPrime,
                                             const linear::Vector& b)
    {
        // Same as residualNorm but specialized for pressure correction
        // Computes ||A p' - b||_2 where A is pressure-correction matrix
        return residualNorm(A, pPrime, b);
    }

    // Turbulence field residual: L2 norm of field difference ||φ_new - φ_old||
    // Used to track convergence of k, ε fields
    static linear::Scalar fieldResidualL2(const core::Field& field_new,
                                        const core::Field& field_old)
    {
        if (field_new.size() != field_old.size()) {
            throw std::runtime_error("Field sizes do not match in fieldResidualL2");
        }

        linear::Scalar sum = 0.0;
        for (int idx = 0; idx < field_new.size(); ++idx) {
            linear::Scalar diff = field_new[idx] - field_old[idx];
            sum += diff * diff;
        }

        return (field_new.size() > 0) 
            ? std::sqrt(sum / field_new.size()) 
            : 0.0;
    }

    // Convenience methods for turbulence fields
    static linear::Scalar turbulentKineticEnergyResidualL2(
        const core::Field& k_new,
        const core::Field& k_old)
    {
        return fieldResidualL2(k_new, k_old);
    }

    static linear::Scalar turbulentDissipationResidualL2(
        const core::Field& epsilon_new,
        const core::Field& epsilon_old)
    {
        return fieldResidualL2(epsilon_new, epsilon_old);
    }

    // Production-to-dissipation ratio for k-equation source balance check
    // Returns average Pk / (ρ ε) over interior domain
    // Value > 1 indicates production dominated; < 1 indicates dissipation dominated
    static linear::Scalar productionDissipationRatio(
        const core::Field& production,  // Pk
        const core::Field& k,
        const core::Field& epsilon,
        linear::Scalar rho = 1.0)
    {
        const int Nx = k.size();  // Assuming scalar field indexing
        linear::Scalar sum = 0.0;
        int count = 0;

        for (int idx = 0; idx < Nx; ++idx) {
            // Avoid division by zero
            linear::Scalar eps_safe = std::max(epsilon[idx], 1e-12);
            linear::Scalar diss = rho * eps_safe;
            linear::Scalar prod = production[idx];

            if (diss > 1e-16) {
                sum += prod / diss;
                ++count;
            }
        }

        return (count > 0) ? sum / count : 0.0;
    }

    // Turbulent viscosity ratio: νt/ν average and max
    // Returns {mean_ratio, max_ratio}
    static std::pair<linear::Scalar, linear::Scalar> turbulentViscosityRatio(
        const core::Field& nu_t,
        linear::Scalar nu)
    {
        if (nu < 1e-16) {
            return {0.0, 0.0};
        }

        linear::Scalar sum = 0.0;
        linear::Scalar max_ratio = 0.0;

        for (int idx = 0; idx < nu_t.size(); ++idx) {
            linear::Scalar ratio = nu_t[idx] / nu;
            sum += ratio;
            max_ratio = std::max(max_ratio, ratio);
        }

        linear::Scalar mean_ratio = (nu_t.size() > 0) 
            ? sum / nu_t.size() 
            : 0.0;

        return {mean_ratio, max_ratio};
    }

    // Wall y+ diagnostics: compute max and mean y+ for wall-adjacent cells
    // Returns {max_y_plus, mean_y_plus}
    static std::pair<linear::Scalar, linear::Scalar> wallYPlus(
        const core::UField& u,
        const core::VField& v,
        const core::Mesh& mesh,
        linear::Scalar nu,
        linear::Scalar rho = 1.0)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();
        
        linear::Scalar max_yp = 0.0;
        linear::Scalar sum_yp = 0.0;
        int count = 0;
        
        // Convert to cell-centered for convenience
        core::Field u_cc = u.cellCenterView(mesh);
        core::Field v_cc = v.cellCenterView(mesh);
        
        // Bottom and top walls
        for (int i = 1; i < Nx - 1; ++i) {
            // Bottom wall (j=0): first cell at j=1
            {
                linear::Scalar y_p = dy / 2.0;  // Distance to wall
                linear::Scalar u_p = std::abs(u_cc(i, 1));
                linear::Scalar v_p = std::abs(v_cc(i, 1));
                linear::Scalar vel_mag = std::sqrt(u_p * u_p + v_p * v_p);
                linear::Scalar tau_w = nu * vel_mag / y_p;  // Approximation
                linear::Scalar u_tau = std::sqrt(std::abs(tau_w) / rho);
                linear::Scalar yp = y_p * u_tau / nu;
                max_yp = std::max(max_yp, yp);
                sum_yp += yp;
                ++count;
            }
            
            // Top wall (j=Ny-1): first cell at j=Ny-2
            {
                linear::Scalar y_p = dy / 2.0;
                linear::Scalar u_p = std::abs(u_cc(i, Ny - 2));
                linear::Scalar v_p = std::abs(v_cc(i, Ny - 2));
                linear::Scalar vel_mag = std::sqrt(u_p * u_p + v_p * v_p);
                linear::Scalar tau_w = nu * vel_mag / y_p;
                linear::Scalar u_tau = std::sqrt(std::abs(tau_w) / rho);
                linear::Scalar yp = y_p * u_tau / nu;
                max_yp = std::max(max_yp, yp);
                sum_yp += yp;
                ++count;
            }
        }
        
        // Left and right walls
        for (int j = 1; j < Ny - 1; ++j) {
            // Left wall (i=0): first cell at i=1
            {
                linear::Scalar y_p = dx / 2.0;
                linear::Scalar u_p = std::abs(u_cc(1, j));
                linear::Scalar v_p = std::abs(v_cc(1, j));
                linear::Scalar vel_mag = std::sqrt(u_p * u_p + v_p * v_p);
                linear::Scalar tau_w = nu * vel_mag / y_p;
                linear::Scalar u_tau = std::sqrt(std::abs(tau_w) / rho);
                linear::Scalar yp = y_p * u_tau / nu;
                max_yp = std::max(max_yp, yp);
                sum_yp += yp;
                ++count;
            }
            
            // Right wall (i=Nx-1): first cell at i=Nx-2
            {
                linear::Scalar y_p = dx / 2.0;
                linear::Scalar u_p = std::abs(u_cc(Nx - 2, j));
                linear::Scalar v_p = std::abs(v_cc(Nx - 2, j));
                linear::Scalar vel_mag = std::sqrt(u_p * u_p + v_p * v_p);
                linear::Scalar tau_w = nu * vel_mag / y_p;
                linear::Scalar u_tau = std::sqrt(std::abs(tau_w) / rho);
                linear::Scalar yp = y_p * u_tau / nu;
                max_yp = std::max(max_yp, yp);
                sum_yp += yp;
                ++count;
            }
        }
        
        linear::Scalar mean_yp = (count > 0) ? sum_yp / count : 0.0;
        return {max_yp, mean_yp};
    }
};

}

