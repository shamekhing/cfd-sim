#include "../../include/linear/BiCGSTAB.hpp"
#include "../../backend/BackendDispatcher.hpp"
#include <iostream>
#include <stdexcept>
#include <chrono>

namespace linear {

BiCGSTABSolver::BiCGSTABSolver(int maxIter, Scalar tol) 
    : maxIter_(maxIter), tol_(tol) {}

Vector BiCGSTABSolver::solve(const SparseMatrix& A, const Vector& b) {
    try {
        auto start = std::chrono::high_resolution_clock::now();
        
        int iters = 0;
        ::BackendImpl::Vector backend_result = ::BackendImpl::solveBiCGSTAB(A.raw(), b.raw(), maxIter_, tol_, iters);
        
        Vector x;
        x.setRawData(backend_result);
        
        auto end = std::chrono::high_resolution_clock::now();
        solveTime_ = std::chrono::duration<double, std::milli>(end - start).count();
        
        iters_ = iters;
        residual_ = computeResidual(A, x, b);
        converged_ = true;
        
        if (verbosity_ > 0) {
            std::cout << "BiCGSTAB: " << iters_ << " iterations, res=" << residual_ 
                      << ", time=" << solveTime_ << "ms" << std::endl;
        }
        
        return x;
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("BiCGSTAB solver error: ") + e.what());
    }
}

}  // namespace linear
