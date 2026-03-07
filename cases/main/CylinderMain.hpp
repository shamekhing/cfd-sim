#pragma once
/**
 * @file CylinderMain.hpp
 * @brief Flagship vortex-shedding case – cylinder dead-centre in the domain.
 *
 * Domain : [0, Lx] × [0, Ly]
 * Cylinder: centred at (Lx/2, Ly/2), radius r.
 *
 * Boundary conditions (properly imposed by overriding the assembler RHS):
 *   Left  (x=0)    : u = U_inf, v = 0  (uniform inlet)
 *   Right (x=Lx)   : du/dx = 0, v = 0  (zero-gradient outlet)
 *   Top/Bottom      : u = U_inf, v = 0  (free-slip / free-stream)
 *   Cylinder surface: u = v = 0         (no-slip, IBM mask)
 *
 * Physics target (Williamson 1989): Re=100 → St≈0.164, Cd≈1.3
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>

#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"
#include "../../include/core/UField.hpp"
#include "../../include/core/VField.hpp"
#include "../../include/core/Indexing.hpp"
#include "../../include/linear/Vector.hpp"
#include "../../include/linear/BiCGSTAB.hpp"
#include "../../include/assembly/StaggeredMomentumAssembler.hpp"
#include "../../include/assembly/StaggeredPressureCorrectionAssembler.hpp"
#include "../../include/io/PVDWriter.hpp"
#include "../../include/io/VTSReader.hpp"
#include "../../include/config/CaseConfigs.hpp"

namespace cases {

// ─────────────────────────────────────────────────────────────────────────────
// IBM / BC helpers (all in anonymous detail namespace)
// ─────────────────────────────────────────────────────────────────────────────
namespace cyl_main_detail {

/// True if the u-face at staggered position (i*dx, (j+0.5)*dy) is inside the cylinder.
inline bool uFaceInCylinder(int i, int j,
                             const core::Mesh& mesh,
                             double cx, double cy, double r)
{
    double x = i * mesh.dx();
    double y = (j + 0.5) * mesh.dy();
    double dx = x - cx, dy = y - cy;
    return (dx*dx + dy*dy) <= r*r;
}

/// True if the v-face at staggered position ((i+0.5)*dx, j*dy) is inside the cylinder.
inline bool vFaceInCylinder(int i, int j,
                             const core::Mesh& mesh,
                             double cx, double cy, double r)
{
    double x = (i + 0.5) * mesh.dx();
    double y = j * mesh.dy();
    double dx = x - cx, dy = y - cy;
    return (dx*dx + dy*dy) <= r*r;
}

/// Zero all velocity faces that fall inside the cylinder (IBM no-slip).
inline void maskCylinder(core::UField& u, core::VField& v,
                          const core::Mesh& mesh,
                          double cx, double cy, double r)
{
    const int Nx = mesh.Nx(), Ny = mesh.Ny();
    for (int j = 0; j < Ny; ++j)
        for (int i = 0; i <= Nx; ++i)
            if (uFaceInCylinder(i, j, mesh, cx, cy, r)) u(i, j) = 0.0;
    for (int j = 0; j <= Ny; ++j)
        for (int i = 0; i < Nx; ++i)
            if (vFaceInCylinder(i, j, mesh, cx, cy, r)) v(i, j) = 0.0;
}

/**
 * Override boundary rows in the u-momentum RHS to impose cylinder-flow BCs.
 * The assembler already wrote lid values; we replace them.
 *
 *  Left  (i==0     ): b = U_inf   (inlet)
 *  Right (i==Nx    ): b = u(Nx-1,j) (zero-gradient outlet)
 *  Bottom(j==0     ): b = U_inf   (free-slip)
 *  Top   (j==Ny-1  ): b = U_inf   (free-slip, was lid=1 → same if U_inf=1)
 */
inline void overrideUBC(linear::Vector& b,
                         const core::Mesh& mesh,
                         const core::UField& u,
                         double U_inf)
{
    const int Nx = mesh.Nx(), Ny = mesh.Ny();
    using core::idxU;

    // Left inlet and right outlet (full column of u-faces)
    for (int j = 0; j < Ny; ++j) {
        b(idxU(0,  j, Nx)) = static_cast<linear::Scalar>(U_inf);        // inlet
        b(idxU(Nx, j, Nx)) = u(Nx - 1, j);                              // outlet zero-grad
    }
    // Bottom and top slip rows
    for (int i = 0; i <= Nx; ++i) {
        b(idxU(i, 0,      Nx)) = static_cast<linear::Scalar>(U_inf);    // bottom
        b(idxU(i, Ny - 1, Nx)) = static_cast<linear::Scalar>(U_inf);    // top
    }
}

/**
 * Override boundary rows in the v-momentum RHS.
 * The assembler sets them all to 0 already (lid BC for v is 0); for a
 * cylinder-flow case v=0 on all boundaries so no change is needed —
 * this function is a no-op but kept for clarity / future extension.
 */
inline void overrideVBC(linear::Vector& /*b*/, const core::Mesh& /*mesh*/) {}

/**
 * Compute drag (Cd) and lift (Cl) coefficients by integrating pressure and
 * viscous stress over cells in a thin annular band around the cylinder.
 */
inline std::pair<double, double> computeForceCoefficients(
    const core::Mesh& mesh,
    const core::Field& p,
    const core::UField& u,
    const core::VField& v,
    double cx, double cy, double r,
    double U_inf, double nu)
{
    const int Nx = mesh.Nx(), Ny = mesh.Ny();
    const double dx = mesh.dx(), dy = mesh.dy();
    const double D  = 2.0 * r;

    // Surface cells: pressure cells whose centre is in the one-cell-wide
    // annular band just outside the cylinder (inner = r, outer = r + dx).
    //
    // Normalization: each surface cell is assigned an equal share of the
    // cylinder circumference (π·D / N_surface), so the total arc length=π·D
    // regardless of grid resolution. This keeps Cd bounded on coarse grids
    // where D/dx can be as small as 3-5 cells.
    const double inner = r;
    const double outer = r + std::max(dx, dy);

    // Pass 1: count surface cells
    int Ns = 0;
    for (int j = 1; j < Ny - 1; ++j)
        for (int i = 1; i < Nx - 1; ++i) {
            double ddx = mesh.x(i) - cx, ddy = mesh.y(j) - cy;
            double dist = std::sqrt(ddx*ddx + ddy*ddy);
            if (dist > inner && dist <= outer) ++Ns;
        }
    if (Ns == 0) return {0.0, 0.0};

    const double dS = M_PI * D / static_cast<double>(Ns);   // arc per cell

    // Pass 2: integrate pressure × outward normal × dS
    double Fx = 0.0, Fy = 0.0;
    for (int j = 1; j < Ny - 1; ++j) {
        for (int i = 1; i < Nx - 1; ++i) {
            double ddx = mesh.x(i) - cx, ddy = mesh.y(j) - cy;
            double dist = std::sqrt(ddx*ddx + ddy*ddy);
            if (dist <= inner || dist > outer) continue;

            double nx = ddx / dist, ny_n = ddy / dist;

            // Pressure contribution
            Fx -= p(i, j) * nx   * dS;
            Fy -= p(i, j) * ny_n * dS;

            // Viscous (first-order approximate)
            double dudy = (u(i + 1, j) - u(i, j)) / dx;
            double dvdx = (v(i, j + 1) - v(i, j)) / dy;
            Fx += nu * dudy * ny_n * dS;
            Fy += nu * dvdx * nx   * dS;
        }
    }

    double qD = 0.5 * U_inf * U_inf * D;   // reference force per unit depth
    return { Fx / qD, Fy / qD };
}

/**
 * Estimate Strouhal number from lift-coefficient history via zero-crossing
 * period detection.
 */
inline double strouhalFromCl(const std::vector<double>& cl_hist,
                               double dt, double U_inf, double D,
                               size_t startIndex = 0)
{
    if (cl_hist.size() < 20 || startIndex >= cl_hist.size()) return 0.0;
    const size_t n = cl_hist.size();
    const size_t count = n - startIndex;
    if (count < 20) return 0.0;

    double meanCl = 0.0;
    for (size_t k = startIndex; k < n; ++k) meanCl += cl_hist[k];
    meanCl /= static_cast<double>(count);

    std::vector<int> crossings;
    for (size_t k = std::max<size_t>(startIndex + 1, 1); k < n; ++k) {
        const double prev = cl_hist[k - 1] - meanCl;
        const double curr = cl_hist[k] - meanCl;
        if (prev < 0.0 && curr >= 0.0)
            crossings.push_back(static_cast<int>(k));
    }
    if (crossings.size() < 2) return 0.0;
    double T = static_cast<double>(crossings.back() - crossings[crossings.size() - 2]) * dt;
    return D / (U_inf * T);
}

}  // namespace cyl_main_detail


// ─────────────────────────────────────────────────────────────────────────────
// runCylinderMain — public entry point
// ─────────────────────────────────────────────────────────────────────────────
inline void runCylinderMain(
    const config::CylinderFlowConfig& cfg,
    const std::string& outputDir = "build/output",
    const std::string& name      = "cylinder")
{
    using namespace core;
    namespace det = cyl_main_detail;

    const int    Nx      = cfg.Nx;
    const int    Ny      = cfg.Ny;
    const double Lx      = cfg.Lx;
    const double Ly      = cfg.Ly;
    const double U_inf   = cfg.U_inf;
    const double r       = cfg.radius;
    const double nu      = cfg.nu();
    const double dt      = cfg.dt;
    const int    nSteps  = cfg.numTimeSteps;
    const double D       = 2.0 * r;

    // ── Cylinder dead-centre ──────────────────────────────────────────────────
    const double cx = Lx / 2.0;
    const double cy = Ly / 2.0;

    // ── Output directories ───────────────────────────────────────────────────
    const std::string runDir = outputDir + "/" + name;
    const std::string vtsDir = runDir + "/vts";
    const std::string csvDir = runDir + "/csv";
    std::filesystem::create_directories(vtsDir);
    std::filesystem::create_directories(csvDir);

    if (cfg.verbosity >= 1) {
        std::cout << "\n╔══════════════════════════════════════════╗\n";
        std::cout << "║  Cylinder Main  –  Flagship Case         ║\n";
        std::cout << "╚══════════════════════════════════════════╝\n";
        std::cout << "  Domain   : [0," << Lx << "] × [0," << Ly << "]\n";
        std::cout << "  Cylinder : centre (" << cx << "," << cy
                  << ")  D=" << D << "  r=" << r << "\n";
        std::cout << "  Re=" << cfg.Re << "  ν=" << nu
                  << "  U_inf=" << U_inf << "\n";
        std::cout << "  Grid     : " << Nx << " × " << Ny << "\n";
        std::cout << "  dt=" << dt << "  steps=" << nSteps << "\n";
        std::cout << "  SIMPLE   : " << cfg.itersPerTimeStep
                  << " inner iters/step\n";
        std::cout << "  Ref      : St≈0.164, Cd≈1.3  "
                  << "(Williamson 1989, Re=100)\n\n";
    }

    // ── Mesh and fields ──────────────────────────────────────────────────────
    Mesh   mesh(Nx, Ny, Lx, Ly);
    UField u(mesh), u_old(mesh);
    VField v(mesh), v_old(mesh);
    Field  p(mesh);

    u.setZero(); v.setZero(); p.setZero();

    // ── PVD writer with RESUME support ───────────────────────────────────────
    io::PVDWriter pvd(vtsDir, name, cfg.verbosity, runDir);
    const int startFrame     = pvd.getStartFrame();
    const std::string lastVts = pvd.getLastFrameFilename();

    // Resume: read the last written cell-centre fields back into staggered arrays
    if (startFrame > 0 && !lastVts.empty()) {
        Field uc(mesh), vc(mesh), pc(mesh);
        if (io::VTSReader::read(lastVts, mesh, uc, vc, pc)) {
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i) {
                    u(i, j) = uc(i, j);
                    v(i, j) = vc(i, j);
                    p(i, j) = pc(i, j);
                }
            if (cfg.verbosity >= 1)
                std::cout << "  Resumed from frame " << startFrame << "\n";
        } else {
            std::cerr << "  Warning: could not read " << lastVts
                      << "; starting from t=0\n";
        }
    }

    // Initial condition: free-stream everywhere outside cylinder
    if (startFrame == 0) {
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i <= Nx; ++i)
                if (!det::uFaceInCylinder(i, j, mesh, cx, cy, r))
                    u(i, j) = static_cast<linear::Scalar>(U_inf);

        // One-shot antisymmetric v impulse in the near-wake:
        // seeds the shear layer asymmetry at t=0 without waiting for
        // the inlet perturbation to convect to the cylinder (~10 dt).
        const double amp = 0.05 * U_inf;   // 5% — decays quickly, just seeds the mode
        const double mdx = mesh.dx(), mdy = mesh.dy();
        for (int j = 1; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double xi = (i + 0.5) * mdx;
                double yj = j * mdy;
                double ddx = xi - cx, ddy = yj - cy;
                double dist = std::sqrt(ddx*ddx + ddy*ddy);
                // Only cells immediately downstream of the cylinder
                if (ddx > 0.0 && dist < 4.0 * r) {
                    // Antisymmetric: positive above centreline, negative below
                    double sign = (yj > cy) ? 1.0 : -1.0;
                    v(i, j) += static_cast<linear::Scalar>(sign * amp
                        * std::exp(-dist / (2.0 * r)));
                }
            }
        }
        det::maskCylinder(u, v, mesh, cx, cy, r);   // clean up any masked faces
    }

    // ── CSV output ───────────────────────────────────────────────────────────
    std::ofstream fForces(csvDir + "/forces.csv");
    fForces << "step,time,Cd,Cl\n";

    std::ofstream fSt(csvDir + "/strouhal.csv");
    fSt << "step,time,f_shed,St\n";

    // ── Solver ───────────────────────────────────────────────────────────────
    linear::BiCGSTABSolver lsolver(cfg.maxLinearIter,
                                    static_cast<linear::Scalar>(cfg.linearTol));

    std::vector<double> cl_hist;
    cl_hist.reserve(nSteps);
    double lastSt = 0.0;

    const double alphaP  = static_cast<double>(cfg.alpha_p);
    const double dx      = mesh.dx();
    const double dy      = mesh.dy();

    auto t0 = std::chrono::high_resolution_clock::now();

    // ══════════════════════════════════════════════════════════════════════════
    // Time-stepping loop
    // ══════════════════════════════════════════════════════════════════════════
    for (int step = startFrame; step < nSteps; ++step) {

        const double time = (step + 1) * dt;
        const double rampFactor = (time < 2.0)
            ? 0.5 * (1.0 - std::cos(M_PI * time / 2.0))
            : 1.0;
        const double U_step = U_inf * rampFactor;
        u_old = u;
        v_old = v;

        // ── SIMPLE inner iterations ──────────────────────────────────────────
        // For *transient* SIMPLE only 2-5 iterations per step are correct.
        // Over-converging (e.g. 200 iters) acts like a pseudo-steady solver
        // and suppresses physical unsteadiness.
        for (int it = 0; it < cfg.itersPerTimeStep; ++it) {

            // ─ U momentum ────────────────────────────────────────────────────
            auto sysU = assembly::StaggeredMomentumAssembler::assembleU(
                mesh, u, v, p,
                static_cast<linear::Scalar>(nu),
                static_cast<linear::Scalar>(dt),
                &u_old);

            det::overrideUBC(sysU.b, mesh, u, U_step);

            linear::Vector uStar = lsolver.solve(sysU.A, sysU.b);

            // ─ V momentum ────────────────────────────────────────────────────
            auto sysV = assembly::StaggeredMomentumAssembler::assembleV(
                mesh, u, v, p,
                static_cast<linear::Scalar>(nu),
                static_cast<linear::Scalar>(dt),
                &v_old);

            det::overrideVBC(sysV.b, mesh);   // no-op: already correct

            linear::Vector vStar = lsolver.solve(sysV.A, sysV.b);

            // ─ Write u*, v* back to staggered fields ──────────────────────────
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i <= Nx; ++i)
                    u(i, j) = uStar(idxU(i, j, Nx));

            for (int j = 0; j <= Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    v(i, j) = vStar(idxV(i, j, Nx));

            // ─ Cylinder IBM mask ──────────────────────────────────────────────
            det::maskCylinder(u, v, mesh, cx, cy, r);

            // ─ Pressure correction ────────────────────────────────────────────
            auto sysP = assembly::StaggeredPressureCorrectionAssembler::assemble(
                mesh, u, v, sysU.faceDiagU, sysV.faceDiagV);

            linear::Vector pCorr = lsolver.solve(sysP.A, sysP.b);

            // ─ Velocity correction ────────────────────────────────────────────
            // u-face (i,j): east face of cell (i-1,j), interior i in [1..Nx-1]
            for (int j = 0; j < Ny; ++j) {
                for (int i = 1; i < Nx; ++i) {
                    int row  = idxU(i, j, Nx);
                    double ae = static_cast<double>(sysU.faceDiagU(row));
                    if (std::abs(ae) < 1e-30) continue;
                    double de = dy / ae;
                    u(i, j) += static_cast<linear::Scalar>(
                        de * (pCorr(id(i - 1, j, Nx)) - pCorr(id(i, j, Nx))));
                }
            }
            // v-face (i,j): north face of cell (i,j-1), interior j in [1..Ny-1]
            for (int j = 1; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    int row  = idxV(i, j, Nx);
                    double an = static_cast<double>(sysV.faceDiagV(row));
                    if (std::abs(an) < 1e-30) continue;
                    double dn = dx / an;
                    v(i, j) += static_cast<linear::Scalar>(
                        dn * (pCorr(id(i, j - 1, Nx)) - pCorr(id(i, j, Nx))));
                }
            }

            // ─ Pressure update ────────────────────────────────────────────────
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    p(i, j) += static_cast<linear::Scalar>(
                        alphaP * static_cast<double>(pCorr(id(i, j, Nx))));

            // ─ Re-apply BCs and cylinder mask ─────────────────────────────────
            // Benchmark-style uniform inlet with smooth startup ramp.
            for (int j = 0; j < Ny; ++j) u(0, j) = static_cast<linear::Scalar>(U_step);
            for (int j = 0; j < Ny; ++j)  u(Nx, j) = u(Nx - 1, j);
            for (int i = 0; i <= Nx; ++i) u(i, 0)      = static_cast<linear::Scalar>(U_step);
            for (int i = 0; i <= Nx; ++i) u(i, Ny - 1) = static_cast<linear::Scalar>(U_step);
            det::maskCylinder(u, v, mesh, cx, cy, r);

        }  // SIMPLE iterations

        // ── Cell-centre views for output ──────────────────────────────────────
        Field u_cc = u.cellCenterView(mesh);
        Field v_cc = v.cellCenterView(mesh);

        // ── Forces ───────────────────────────────────────────────────────────
        auto [Cd, Cl] = det::computeForceCoefficients(
            mesh, p, u, v, cx, cy, r, U_inf, nu);
        cl_hist.push_back(Cl);
        fForces << step + 1 << "," << std::fixed << std::setprecision(6)
                << time << "," << Cd << "," << Cl << "\n";

        // ── Strouhal ─────────────────────────────────────────────────────────
        const size_t n = cl_hist.size();
        const size_t windowSamples = std::max<size_t>(20, static_cast<size_t>(100.0 / dt));
        const size_t startIndex = (n > windowSamples) ? (n - windowSamples) : (n / 2);
        double St = det::strouhalFromCl(cl_hist, dt, U_inf, D, startIndex);
        if (St > 0.05 && St < 0.50) {
            if (cfg.verbosity >= 1 && std::abs(St - lastSt) > 0.005)
                std::cout << "\n  ★  Strouhal: St=" << std::fixed
                          << std::setprecision(4) << St << "  (ref 0.164)\n";
            fSt << step + 1 << "," << time << ","
                << St * U_inf / D << "," << St << "\n";
        }
        lastSt = St;

        // ── Progress bar ─────────────────────────────────────────────────────
        {
            auto now = std::chrono::high_resolution_clock::now();
            double el = std::chrono::duration<double>(now - t0).count();
            double eta = el / static_cast<double>(step - startFrame + 1)
                           * static_cast<double>(nSteps - step - 1);
            int pct = (step + 1) * 100 / nSteps;
            std::cout << "\r  ["
                      << std::string(pct / 5, '#')
                      << std::string(20 - pct / 5, '-')
                      << "] " << std::setw(3) << pct
                      << "% step " << step + 1 << "/" << nSteps
                      << " T=" << std::fixed << std::setprecision(1) << el
                      << "s ETA=" << eta << "s"
                      << " Cd=" << std::setprecision(3) << Cd
                      << " Cl=" << Cl
                      << std::flush;
        }

        // ── Animation frame ──────────────────────────────────────────────────
        if (cfg.useAnimation)
            pvd.writeTimestep(step, time, mesh, u_cc, v_cc, p);

    }  // time loop

    // ── Finalise ─────────────────────────────────────────────────────────────
    std::cout << "\n";
    fForces.close();
    fSt.close();

    if (cfg.useAnimation) {
        pvd.finalize();
        pvd.writePVD();
    }

    if (cfg.verbosity >= 1) {
        std::cout << "\n======================================\n";
        std::cout << "  Cylinder Main  —  complete\n";
        std::cout << "  Final St   : " << lastSt << "\n";
        std::cout << "  Ref St     : 0.164  (Williamson 1989)\n";
        std::cout << "  Output     : " << runDir << "\n";
        std::cout << "======================================\n";
    }
}

}  // namespace cases
