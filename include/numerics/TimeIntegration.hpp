#pragma once

#include "../linear/Vector.hpp"

namespace numerics {

/**
 * @brief Time integration schemes for transient problems
 * 
 * Currently empty - supports future expansion to:
 * - Explicit Euler (1st order, conditionally stable)
 * - Implicit Euler (1st order, A-stable)
 * - Crank-Nicolson (2nd order, A-stable)
 * - BDF2, BDF3 (higher order, stiff)
 * - RK2, RK3, RK4 (explicit Runge-Kutta methods)
 */
class TimeIntegration {
public:
    // Future implementations
};

}  // namespace numerics
