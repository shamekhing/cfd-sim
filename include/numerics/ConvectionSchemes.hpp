#pragma once

#include "../linear/Vector.hpp"
#include <cmath>
#include <algorithm>

namespace numerics {

enum class ConvectionScheme {
    Central,   // 2nd-order, non-monotonic
    Upwind,    // 1st-order, monotonic
    QUICK,     // 3rd-order
    MUSCL,     // 2nd-order TVD
    VanLeer    // 2nd-order TVD, smooth limiter
};

/**
 * @brief 5 convection interpolation schemes for finite volume method
 * 
 * Interpolates scalar φ to face for flux computation:
 * - Central: Unbounded, 2nd-order, prone to oscillations
 * - Upwind: 1st-order, bounded, stable but dissipative
 * - QUICK: 3rd-order, bounded with limiters
 * - MUSCL: 2nd-order TVD with monotone limiter
 * - VanLeer: 2nd-order TVD with smooth limiter
 */
class ConvectionSchemes {
public:
    /**
     * @brief Interpolate scalar φ to face using specified scheme
     * 
     * @param phi_far_upstream   φ value 2 cells upstream
     * @param phi_upstream       φ value 1 cell upstream
     * @param phi_center         φ value at center (current cell)
     * @param phi_downstream     φ value 1 cell downstream
     * @param velocity           Local velocity (flow direction)
     * @param grid_spacing       Grid spacing (dx or dy)
     * @param scheme             Discretization scheme
     * @return                   Interpolated value at face
     */
    static linear::Scalar interpolateFace(
        linear::Scalar phi_far_upstream,
        linear::Scalar phi_upstream,
        linear::Scalar phi_center,
        linear::Scalar phi_downstream,
        linear::Scalar velocity,
        linear::Scalar grid_spacing,
        ConvectionScheme scheme);

    // Flux limiters for TVD schemes
    static linear::Scalar vanLeerLimiter(linear::Scalar r);
    static linear::Scalar musclLimiter(linear::Scalar r);
    static linear::Scalar minmodLimiter(linear::Scalar r);
    static linear::Scalar superbeeLimiter(linear::Scalar r);

private:
    static linear::Scalar computeSlope(linear::Scalar phi_up, linear::Scalar phi_center, linear::Scalar phi_down);
};

}  // namespace numerics
