#pragma once

#include "Solver.hpp"

namespace linear {

class LUSolver : public Solver {
private:
    int iters_ = 1;  // Direct solvers take 1 "iteration"
    Scalar residual_ = 0;
    bool converged_ = true;
    int verbosity_ = 0;

public:
    LUSolver() = default;

    Vector solve(const SparseMatrix& A, const Vector& b) override;

    int iterations() const override { return iters_; }
    Scalar residual() const override { return residual_; }
    bool converged() const override { return converged_; }

    void setMaxIterations(int maxIter) override {}  // ignored for direct solver
    void setTolerance(Scalar tol) override {}        // ignored for direct solver
    void setVerbosity(int verbosity) override { verbosity_ = verbosity; }

    std::string name() const override { return "LU"; }
};

}  // namespace linear
