#pragma once

#include "../include/linear/BiCGSTAB.hpp"
#include "../include/linear/LU.hpp"
#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "../include/core/UField.hpp"
#include "../include/core/VField.hpp"
#include "../include/core/Indexing.hpp"
#include "../include/assembly/StaggeredMomentumAssembler.hpp"
#include "../include/assembly/StaggeredPressureCorrectionAssembler.hpp"
#include "../include/physics/TurbulenceKEpsilon.hpp"
#include <iostream>
#include <stdexcept>
#include <string>
#include <chrono>
#include <iomanip>

namespace solver {

class StaggeredSIMPLE {
public:
    static void applyVelocityBC(core::UField& u, core::VField& v, int Nx, int Ny)
    {
        for (int i = 0; i <= Nx; ++i) {
            u(i, 0) = 0.0;
            u(i, Ny - 1) = 1.0;
        }
        for (int j = 0; j < Ny; ++j) {
            u(0, j) = 0.0;
            u(Nx, j) = 0.0;
        }
        for (int i = 0; i < Nx; ++i) {
            v(i, 0) = 0.0;
            v(i, Ny) = 0.0;
        }
        for (int j = 0; j <= Ny; ++j) {
            v(0, j) = 0.0;
            v(Nx - 1, j) = 0.0;
        }
    }

    static void solveSteady(
        core::Mesh& mesh,
        core::UField& u,
        core::VField& v,
        core::Field& p,
        linear::Scalar nu,
        int maxIter,
        linear::Scalar alpha_u = 0.7,
        linear::Scalar alpha_v = 0.7,
        linear::Scalar alpha_p = 0.3,
        int verbosity = 1,
        physics::TurbulenceKEpsilon* turbulence = nullptr,
        linear::Scalar source_u = 0.0,
        linear::Scalar source_v = 0.0)
    {
        using namespace core;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        // Quasi-transient SIMPLE: use fictitious dt for diagonal stability
        // Courant-like parameter: dt = 0.1 * dx to add 10% diagonal term without over-damping
        const linear::Scalar dt_pseudo = static_cast<linear::Scalar>(0.1) * dx;

        applyVelocityBC(u, v, Nx, Ny);

        // Start timing
        auto start_time = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < maxIter; ++iter) {
            // Calculate elapsed and estimated remaining time
            auto current_time = std::chrono::high_resolution_clock::now();
            double elapsed_sec = std::chrono::duration<double>(current_time - start_time).count();
            double avg_iter_time = (iter > 0) ? elapsed_sec / iter : 0.0;
            int remaining_iters = maxIter - iter;
            double estimated_remaining = avg_iter_time * remaining_iters;

            // Progress output on single line with carriage return
            if (verbosity > 0) {
                int progress = (iter * 100) / (maxIter - 1);
                std::cout << "\r  SIMPLE iter " << iter << "/" << maxIter 
                          << " [" << std::string(progress / 5, '#') 
                          << std::string(20 - progress / 5, '-') << "] " 
                          << progress << "% | " 
                          << std::fixed << std::setprecision(1) 
                          << "T:" << elapsed_sec << "s | ETA:" << estimated_remaining << "s"
                          << std::flush;
                if (iter == maxIter - 1) std::cout << std::endl;
            }

            auto sysU = assembly::StaggeredMomentumAssembler::assembleU(mesh, u, v, p, nu, dt_pseudo, &u, source_u);
            auto sysV = assembly::StaggeredMomentumAssembler::assembleV(mesh, u, v, p, nu, dt_pseudo, &v, source_v);

            linear::Vector solU;
            linear::Vector solV;
            try {
                solU = linear::LUSolver().solve(sysU.A, sysU.b);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("SIMPLE steady: U-momentum solve failed at iter ") + std::to_string(iter) + ": " + e.what());
            }
            try {
                solV = linear::LUSolver().solve(sysV.A, sysV.b);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("SIMPLE steady: V-momentum solve failed at iter ") + std::to_string(iter) + ": " + e.what());
            }

            core::UField ustar(mesh);
            core::VField vstar(mesh);
            ustar.raw() = solU;
            vstar.raw() = solV;
            applyVelocityBC(ustar, vstar, Nx, Ny);

            linear::Vector faceDiagU = sysU.faceDiagU;
            linear::Vector faceDiagV = sysV.faceDiagV;
            auto pSys = assembly::StaggeredPressureCorrectionAssembler::assemble(
                mesh, ustar, vstar, faceDiagU, faceDiagV);

            core::Field pCorr(mesh);
            linear::Vector solP;
            try {
                solP = linear::LUSolver().solve(pSys.A, pSys.b);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("SIMPLE steady: pressure-correction solve failed at iter ") + std::to_string(iter) + ": " + e.what());
            }
            pCorr.raw() = solP;

            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i <= Nx; ++i) {
                    if (core::isUBoundary(i, j, Nx, Ny)) continue;
                    linear::Scalar a_e = faceDiagU(core::idxU(i, j, Nx));
                    if (a_e < static_cast<linear::Scalar>(1e-30)) continue;
                    linear::Scalar d_e = dy / a_e;
                    linear::Scalar pP = pCorr(i - 1, j);
                    linear::Scalar pE = pCorr(i, j);
                    linear::Scalar u_corr = d_e * (pP - pE);
                    linear::Scalar u_new = ustar(i, j) + u_corr;
                    u(i, j) = (1.0 - alpha_u) * u(i, j) + alpha_u * u_new;
                }
            }
            for (int j = 0; j <= Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    if (core::isVBoundary(i, j, Nx, Ny)) continue;
                    linear::Scalar a_n = faceDiagV(core::idxV(i, j, Nx));
                    if (a_n < static_cast<linear::Scalar>(1e-30)) continue;
                    linear::Scalar d_n = dx / a_n;
                    linear::Scalar pP = pCorr(i, j - 1);
                    linear::Scalar pN = pCorr(i, j);
                    linear::Scalar v_corr = d_n * (pP - pN);
                    linear::Scalar v_new = vstar(i, j) + v_corr;
                    v(i, j) = (1.0 - alpha_v) * v(i, j) + alpha_v * v_new;
                }
            }
            p.raw() += alpha_p * pCorr.raw();
            applyVelocityBC(u, v, Nx, Ny);

            // Solve turbulence equations if enabled
            if (turbulence != nullptr) {
                turbulence->solve(u, v, mesh, nu, dt_pseudo);
            }
        }

        // Print final timing summary
        auto end_time = std::chrono::high_resolution_clock::now();
        double total_sec = std::chrono::duration<double>(end_time - start_time).count();
        double avg_iter_time = total_sec / maxIter;
        if (verbosity > 0) {
            std::cout << "  Steady SIMPLE: Total time = " << std::fixed << std::setprecision(2) 
                      << total_sec << "s, Avg per iter = " << avg_iter_time << "s" << std::endl;
        }
    }

    static void solveTransient(
        core::Mesh& mesh,
        core::UField& u,
        core::VField& v,
        core::Field& p,
        const core::UField& u_old,
        const core::VField& v_old,
        linear::Scalar nu,
        linear::Scalar dt,
        int maxIter,
        linear::Scalar alpha_u = 0.7,
        linear::Scalar alpha_v = 0.7,
        linear::Scalar alpha_p = 0.3,
        int verbosity = 1,
        physics::TurbulenceKEpsilon* turbulence = nullptr)
    {
        using namespace core;

        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const linear::Scalar dx = mesh.dx();
        const linear::Scalar dy = mesh.dy();

        applyVelocityBC(u, v, Nx, Ny);

        // Start timing
        auto start_time = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < maxIter; ++iter) {
            // Calculate elapsed and estimated remaining time
            auto current_time = std::chrono::high_resolution_clock::now();
            double elapsed_sec = std::chrono::duration<double>(current_time - start_time).count();
            double avg_iter_time = (iter > 0) ? elapsed_sec / iter : 0.0;
            int remaining_iters = maxIter - iter;
            double estimated_remaining = avg_iter_time * remaining_iters;

            // Progress output on single line (only if verbosity > 0)
            if (verbosity > 0) {
                int progress = (iter * 100) / (maxIter - 1);
                std::cout << "\r  SIMPLE iter " << iter << "/" << maxIter 
                          << " [" << std::string(progress / 5, '#') 
                          << std::string(20 - progress / 5, '-') << "] " 
                          << progress << "% | "
                          << std::fixed << std::setprecision(1)
                          << "T:" << elapsed_sec << "s | ETA:" << estimated_remaining << "s"
                          << std::flush;
                if (iter == maxIter - 1) std::cout << std::endl;
            }

            auto sysU = assembly::StaggeredMomentumAssembler::assembleU(
                mesh, u, v, p, nu, dt, &u_old);
            auto sysV = assembly::StaggeredMomentumAssembler::assembleV(
                mesh, u, v, p, nu, dt, &v_old);

            linear::Vector solU;
            linear::Vector solV;
            try {
                solU = linear::LUSolver().solve(sysU.A, sysU.b);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("SIMPLE transient: U-momentum solve failed at iter ") + std::to_string(iter) + ": " + e.what());
            }
            try {
                solV = linear::LUSolver().solve(sysV.A, sysV.b);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("SIMPLE transient: V-momentum solve failed at iter ") + std::to_string(iter) + ": " + e.what());
            }

            core::UField ustar(mesh);
            core::VField vstar(mesh);
            ustar.raw() = solU;
            vstar.raw() = solV;
            applyVelocityBC(ustar, vstar, Nx, Ny);

            linear::Vector faceDiagU = sysU.faceDiagU;
            linear::Vector faceDiagV = sysV.faceDiagV;
            auto pSys = assembly::StaggeredPressureCorrectionAssembler::assemble(
                mesh, ustar, vstar, faceDiagU, faceDiagV);

            core::Field pCorr(mesh);
            linear::Vector solP;
            try {
                solP = linear::LUSolver().solve(pSys.A, pSys.b);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("SIMPLE transient: pressure-correction solve failed at iter ") + std::to_string(iter) + ": " + e.what());
            }
            pCorr.raw() = solP;

            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i <= Nx; ++i) {
                    if (core::isUBoundary(i, j, Nx, Ny)) continue;
                    linear::Scalar a_e = faceDiagU(core::idxU(i, j, Nx));
                    if (a_e < static_cast<linear::Scalar>(1e-30)) continue;
                    linear::Scalar d_e = dy / a_e;
                    linear::Scalar pP = pCorr(i - 1, j);
                    linear::Scalar pE = pCorr(i, j);
                    linear::Scalar u_corr = d_e * (pP - pE);
                    linear::Scalar u_new = ustar(i, j) + u_corr;
                    u(i, j) = (1.0 - alpha_u) * u(i, j) + alpha_u * u_new;
                }
            }
            for (int j = 0; j <= Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    if (core::isVBoundary(i, j, Nx, Ny)) continue;
                    linear::Scalar a_n = faceDiagV(core::idxV(i, j, Nx));
                    if (a_n < static_cast<linear::Scalar>(1e-30)) continue;
                    linear::Scalar d_n = dx / a_n;
                    linear::Scalar pP = pCorr(i, j - 1);
                    linear::Scalar pN = pCorr(i, j);
                    linear::Scalar v_corr = d_n * (pP - pN);
                    linear::Scalar v_new = vstar(i, j) + v_corr;
                    v(i, j) = (1.0 - alpha_v) * v(i, j) + alpha_v * v_new;
                }
            }
            p.raw() += alpha_p * pCorr.raw();
            applyVelocityBC(u, v, Nx, Ny);

            // Solve turbulence equations if enabled
            if (turbulence != nullptr) {
                turbulence->solve(u, v, mesh, nu, dt);
            }
        }

        // Print final timing summary
        auto end_time = std::chrono::high_resolution_clock::now();
        double total_sec = std::chrono::duration<double>(end_time - start_time).count();
        double avg_iter_time = total_sec / maxIter;
        if (verbosity > 0) {
            std::cout << "  Transient SIMPLE: Total time = " << std::fixed << std::setprecision(2) 
                      << total_sec << "s, Avg per iter = " << avg_iter_time << "s" << std::endl;
        }
    }

private:
    static void applyPressureBC(core::Field& pCorr, const core::Mesh& mesh)
    {
        int Nx = mesh.Nx();
        int Ny = mesh.Ny();
        for (int j = 0; j < Ny; ++j) {
            pCorr(0, j) = pCorr(1, j);
            pCorr(Nx - 1, j) = pCorr(Nx - 2, j);
        }
        for (int i = 0; i < Nx; ++i) {
            pCorr(i, 0) = pCorr(i, 1);
            pCorr(i, Ny - 1) = pCorr(i, Ny - 2);
        }
    }
};

}
