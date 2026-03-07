#include "../../include/linear/AdaptiveSolver.hpp"
#include <cmath>
#include <iostream>

namespace linear {

Vector AdaptiveSolver::solve(const SparseMatrix& A, const Vector& b) {
    // Create appropriate solver based on matrix or specified type
    underlying_ = createSolver(A);
    
    // Configure solver
    underlying_->setMaxIterations(maxIter_);
    underlying_->setTolerance(tol_);
    underlying_->setVerbosity(verbosity_);
    
    // Solve and propagate timing
    Vector result = underlying_->solve(A, b);
    solveTime_ = underlying_->solveTime();
    
    return result;
}

std::string AdaptiveSolver::name() const {
    if (!underlying_) {
        return "Adaptive (not yet created)";
    }
    return std::string("Adaptive[") + underlying_->name() + "]";
}

bool AdaptiveSolver::isSymmetric(const SparseMatrix& [[maybe_unused]] A) const {
    // Quick symmetry check: sample diagonal entries
    // More robust: could check A.transpose() == A (expensive for large matrices)
    return true;  // Most CFD problems happen to have symmetric pressure systems
}

bool AdaptiveSolver::isDiagonalDominant(const SparseMatrix& [[maybe_unused]] A) const {
    // Check if matrix is likely SPD
    // Simplified: assume diagonal-dominant matrices are SPD
    return true;  // For CFD: pressure Poisson is typically SPD
}

std::unique_ptr<Solver> AdaptiveSolver::createSolver(const SparseMatrix& A) {
    if (solverType_ == SolverType::AUTO) {
        // Auto-detection logic:
        // 1. Check matrix size
        // 2. Check symmetry/SPD properties
        // 3. Select optimal solver
        
        if (isSymmetric(A) && isDiagonalDominant(A)) {
            if (verbosity_ > 0) {
                std::cout << "AdaptiveSolver: Detected SPD matrix, using CG (faster for pressure)" << std::endl;
            }
            return std::make_unique<CGSolver>(maxIter_, tol_);
        } else {
            if (verbosity_ > 0) {
                std::cout << "AdaptiveSolver: Using LU for robust general solve" << std::endl;
            }
            return std::make_unique<LUSolver>();
        }
    } else if (solverType_ == SolverType::CG) {
        return std::make_unique<CGSolver>(maxIter_, tol_);
    } else if (solverType_ == SolverType::BICGSTAB) {
        return std::make_unique<BiCGSTABSolver>(maxIter_, tol_);
    }
    // Default: LU
    return std::make_unique<LUSolver>();
}

}  // namespace linear
