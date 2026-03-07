#include "../../include/linear/LU.hpp"
#include "../../backend/BackendDispatcher.hpp"
#include <iostream>
#include <stdexcept>
#include <chrono>

namespace linear {

Vector LUSolver::solve(const SparseMatrix& A, const Vector& b) {
    try {
        auto start = std::chrono::high_resolution_clock::now();
        
        ::BackendImpl::Vector backend_result = ::BackendImpl::solveLU(A.raw(), b.raw());
        
        Vector x;
        x.setRawData(backend_result);
        
        auto end = std::chrono::high_resolution_clock::now();
        solveTime_ = std::chrono::duration<double, std::milli>(end - start).count();
        
        iters_ = 1;  // Direct solver
        residual_ = computeResidual(A, x, b);
        converged_ = true;
        
        if (verbosity_ > 0) {
            std::cout << "SparseLU: direct solve, res=" << residual_ 
                      << ", time=" << solveTime_ << "ms" << std::endl;
        }
        
        return x;
    } catch (const std::exception& e) {
        converged_ = false;
        residual_ = 1e10;
        iters_ = 0;
        throw std::runtime_error("LU solver failed: " + std::string(e.what()));
    }
}

}  // namespace linear
