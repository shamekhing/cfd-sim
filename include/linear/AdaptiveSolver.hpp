#pragma once

#include "Solver.hpp"
#include "CG.hpp"
#include "BiCGSTAB.hpp"
#include "LU.hpp"
#include <memory>
#include <iostream>

namespace linear {

/**
 * Adaptive Solver - Selects optimal solver based on matrix properties
 * 
 * Strategy:
 * - For pressure Poisson (symmetric, positive-definite): Use CG (fastest)
 * - For momentum equations (general, asymmetric): Use LU (robust, sparse)
 * - Fallback: BiCGSTAB for general matrices
 */
class AdaptiveSolver : public Solver {
public:
    enum class SolverType {
        AUTO,      // Auto-detect based on matrix properties
        CG,        // Conjugate Gradient (symmetric matrices)
        BICGSTAB,  // BiCGSTAB (general asymmetric)
        LU         // Direct LU decomposition (robust)
    };

private:
    std::unique_ptr<Solver> underlying_;
    SolverType solverType_;
    int maxIter_ = 1000;
    Scalar tol_ = 1e-8;
    int verbosity_ = 0;

public:
    AdaptiveSolver(SolverType type = SolverType::AUTO) 
        : solverType_(type), underlying_(nullptr) {}

    Vector solve(const SparseMatrix& A, const Vector& b) override;

    int iterations() const override { 
        return underlying_ ? underlying_->iterations() : 0; 
    }
    Scalar residual() const override { 
        return underlying_ ? underlying_->residual() : 0.0; 
    }
    bool converged() const override { 
        return underlying_ ? underlying_->converged() : false; 
    }

    void setMaxIterations(int maxIter) override { 
        maxIter_ = maxIter;
        if (underlying_) underlying_->setMaxIterations(maxIter);
    }
    
    void setTolerance(Scalar tol) override { 
        tol_ = tol;
        if (underlying_) underlying_->setTolerance(tol);
    }
    
    void setVerbosity(int verbosity) override { 
        verbosity_ = verbosity;
        if (underlying_) underlying_->setVerbosity(verbosity);
    }

    std::string name() const override;

private:
    // Check if matrix is symmetric
    bool isSymmetric(const SparseMatrix& A) const;
    
    // Check if matrix is diagonal-dominant (SPD indicator)
    bool isDiagonalDominant(const SparseMatrix& A) const;
    
    // Create appropriate solver based on matrix type
    std::unique_ptr<Solver> createSolver(const SparseMatrix& A);
};

}  // namespace linear
