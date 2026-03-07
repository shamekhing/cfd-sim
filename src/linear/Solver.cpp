#include "../../include/linear/Solver.hpp"
#include <sstream>
#include <iomanip>

namespace linear {

std::string Solver::info() const {
    std::ostringstream oss;
    oss << name() 
        << ": iters=" << iterations() 
        << ", residual=" << residual() 
        << ", converged=" << (converged() ? "yes" : "no")
        << ", time=" << std::fixed << std::setprecision(3) << solveTime_ << "ms";
    return oss.str();
}

}  // namespace linear
