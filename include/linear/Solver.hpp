#pragma once

#include "../../backend/BackendDispatcher.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include <string>
#include <chrono>

namespace linear {

class Solver {
public:
    virtual ~Solver() = default;

    virtual Vector solve(const SparseMatrix& A, const Vector& b) = 0;

    // Convergence info
    virtual int iterations() const = 0;
    virtual Scalar residual() const = 0;
    virtual bool converged() const = 0;
    
    // Timing info (milliseconds)
    virtual double solveTime() const { return solveTime_; }
    
    // Config
    virtual void setMaxIterations(int maxIter) = 0;
    virtual void setTolerance(Scalar tol) = 0;
    virtual void setVerbosity(int verbosity = 0) = 0;
    
    virtual std::string name() const = 0;
    virtual std::string info() const;

protected:
    double solveTime_ = 0.0;  // Time in milliseconds
    
    // Optimized residual computation (reduces temporary vector creation)
    // Returns ||A*x - b||  with better compiler optimization potential
    inline Scalar computeResidual(const SparseMatrix& A, const Vector& x, const Vector& b) const {
        return (A * x - b).norm();
    }
};

}  // namespace linear
