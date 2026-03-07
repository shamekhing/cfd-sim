#pragma once

#include "Solver.hpp"

namespace linear {

class BiCGSTABSolver : public Solver {
private:
    int maxIter_ = 1000;
    Scalar tol_ = 1e-8;
    int iters_ = 0;
    Scalar residual_ = 0;
    bool converged_ = false;
    int verbosity_ = 0;

public:
    BiCGSTABSolver() = default;
    explicit BiCGSTABSolver(int maxIter, Scalar tol = 1e-8);

    Vector solve(const SparseMatrix& A, const Vector& b) override;

    int iterations() const override { return iters_; }
    Scalar residual() const override { return residual_; }
    bool converged() const override { return converged_; }

    void setMaxIterations(int maxIter) override { maxIter_ = maxIter; }
    void setTolerance(Scalar tol) override { tol_ = tol; }
    void setVerbosity(int verbosity) override { verbosity_ = verbosity; }

    std::string name() const override { return "BiCGSTAB"; }
};

}  // namespace linear
