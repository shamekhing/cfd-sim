#include "../../include/numerics/ConvectionSchemes.hpp"
#include <cmath>
#include <algorithm>

namespace numerics {

linear::Scalar ConvectionSchemes::interpolateFace(
    linear::Scalar phi_far_upstream,
    linear::Scalar phi_upstream,
    linear::Scalar phi_center,
    linear::Scalar phi_downstream,
    linear::Scalar velocity,
    linear::Scalar /*grid_spacing*/,
    ConvectionScheme scheme)
{
    switch (scheme) {
        case ConvectionScheme::Central:
            // Simple average at face
            return 0.5 * (phi_center + phi_downstream);
        
        case ConvectionScheme::Upwind:
            // Choose based on flow direction
            if (velocity >= 0.0) {
                return phi_upstream;
            } else {
                return phi_center;
            }
        
        case ConvectionScheme::QUICK:
            // Quadratic upstream interpolation
            // φ_face = (3/8)φ_up + (3/4)φ_P - (1/8)φ_far_up (for u > 0)
            if (velocity >= 0.0) {
                return (3.0 / 8.0) * phi_upstream +
                       (3.0 / 4.0) * phi_center -
                       (1.0 / 8.0) * phi_far_upstream;
            } else {
                // Downstream side
                return (3.0 / 8.0) * phi_downstream +
                       (3.0 / 4.0) * phi_center -
                       (1.0 / 8.0) * phi_center;
            }
        
        case ConvectionScheme::MUSCL: {
            // Monotonic upwind scheme for conservation laws
            if (velocity >= 0.0) {
                // Upstream side
                linear::Scalar denom = phi_center - phi_upstream;
                if (std::abs(denom) < 1e-16) return phi_upstream;
                linear::Scalar numer = phi_upstream - phi_far_upstream;
                linear::Scalar r = numer / denom;
                linear::Scalar limiter = musclLimiter(r);
                return phi_upstream + 0.5 * limiter * denom;
            } else {
                // Downstream side
                linear::Scalar denom = phi_downstream - phi_center;
                if (std::abs(denom) < 1e-16) return phi_center;
                linear::Scalar numer = phi_center - phi_downstream;
                linear::Scalar r = numer / denom;
                linear::Scalar limiter = musclLimiter(r);
                return phi_center + 0.5 * limiter * denom;
            }
        }
        
        case ConvectionScheme::VanLeer: {
            // Van Leer limiter (smoother than MUSCL)
            if (velocity >= 0.0) {
                linear::Scalar denom = phi_center - phi_upstream;
                if (std::abs(denom) < 1e-16) return phi_upstream;
                linear::Scalar numer = phi_upstream - phi_far_upstream;
                linear::Scalar r = numer / denom;
                linear::Scalar limiter = vanLeerLimiter(r);
                return phi_upstream + 0.5 * limiter * denom;
            } else {
                linear::Scalar denom = phi_downstream - phi_center;
                if (std::abs(denom) < 1e-16) return phi_center;
                linear::Scalar numer = phi_center - phi_downstream;
                linear::Scalar r = numer / denom;
                linear::Scalar limiter = vanLeerLimiter(r);
                return phi_center + 0.5 * limiter * denom;
            }
        }
        
        default:
            return 0.5 * (phi_center + phi_downstream);
    }
}

linear::Scalar ConvectionSchemes::vanLeerLimiter(linear::Scalar r)
{
    // φ(r) = (r + |r|) / (1 + |r|)
    // Ranges [0, 2], smooth, symmetric
    if (r <= 0.0) {
        return 0.0;
    }
    return (r + std::abs(r)) / (1.0 + std::abs(r));
}

linear::Scalar ConvectionSchemes::musclLimiter(linear::Scalar r)
{
    // MUSCL limiter
    // min(max(0, min(2r, (1+r)/2)), 2)
    // Simpler: max(0, min(r, (1+r)/2)) for r >= 0
    if (r <= 0.0) {
        return 0.0;
    }
    return std::max(0.0, std::min(2.0 * r, std::min(1.0, 2.0)));
}

linear::Scalar ConvectionSchemes::minmodLimiter(linear::Scalar r)
{
    // Minmod limiter: min(r, 1) for r > 0
    if (r <= 0.0) {
        return 0.0;
    }
    return std::min(r, 1.0);
}

linear::Scalar ConvectionSchemes::superbeeLimiter(linear::Scalar r)
{
    // Superbee limiter
    if (r <= 0.0) {
        return 0.0;
    }
    if (r >= 1.0) {
        return std::min(r, 2.0);
    }
    return std::max(2.0 * r, std::min(1.0, 2.0));
}

linear::Scalar ConvectionSchemes::computeSlope(linear::Scalar phi_up, linear::Scalar phi_center, linear::Scalar phi_down)
{
    linear::Scalar denom = phi_down - phi_center;
    if (std::abs(denom) < 1e-16) {
        return 0.0;
    }
    linear::Scalar numer = phi_center - phi_up;
    return numer / denom;
}

}  // namespace numerics
