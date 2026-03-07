#pragma once

#include <vector>

/**
 * @brief Reference vortex data for lid-driven cavity validation
 * 
 * Contains primary vortex locations and metrics from literature
 * for different Reynolds numbers.
 */
namespace validation {

/**
 * @brief Vortex center reference data for Re=100
 * 
 * Source: Ghia et al. (1982) "High-Re solutions for incompressible flow..."
 * Journal of Computational Physics
 */
struct VortexRe100 {
    static constexpr double primary_center_x = 0.53;
    static constexpr double primary_center_y = 0.56;
    static constexpr double circulation = 1.89;      // Approximate, normalized
    static constexpr double max_vorticity = 10.5;    // Approximate
    
    // No secondary vortex at Re=100
};

/**
 * @brief Vortex center reference data for Re=400
 * 
 * At Re=400, secondary vortex appears in bottom-left corner
 */
struct VortexRe400 {
    static constexpr double primary_center_x = 0.545;
    static constexpr double primary_center_y = 0.565;
    
    // Secondary vortex in bottom-left corner
    static constexpr double secondary_x = 0.08;
    static constexpr double secondary_y = 0.09;
};

/**
 * @brief Vortex center reference data for Re=1000
 * 
 * At Re=1000, primary vortex shifts and multiple secondary vortices appear.
 * Tertiary vortex may form in top-right corner.
 */
struct VortexRe1000 {
    static constexpr double primary_center_x = 0.545;
    static constexpr double primary_center_y = 0.570;
    
    // Secondary vortex (bottom-left)
    static constexpr double secondary_x = 0.058;
    static constexpr double secondary_y = 0.066;
    
    // Tertiary vortex (top-right, weak)
    static constexpr double tertiary_x = 0.922;
    static constexpr double tertiary_y = 0.918;
};

/**
 * @brief Helper function to get vortex reference for given Re
 * 
 * Returns interpolated or exact reference data for Reynolds number.
 * 
 * @param Re Reynolds number
 * @return Pair of (primary_x, primary_y) centers
 */
inline std::pair<double, double> getVortexReferenceCenter(double Re) {
    if (Re < 50) {
        // Laminar, single centered vortex
        return {0.50, 0.50};
    } else if (Re < 150) {
        // Re=100 region
        return {VortexRe100::primary_center_x, VortexRe100::primary_center_y};
    } else if (Re < 700) {
        // Linear interpolation between Re=100 and Re=400
        // (simplified; literature has more points available)
        double alpha = (Re - 150.0) / (700.0 - 150.0);
        double x = VortexRe100::primary_center_x + 
                   alpha * (VortexRe400::primary_center_x - VortexRe100::primary_center_x);
        double y = VortexRe100::primary_center_y + 
                   alpha * (VortexRe400::primary_center_y - VortexRe100::primary_center_y);
        return {x, y};
    } else {
        // Re >= 700: use Re=1000 data
        return {VortexRe1000::primary_center_x, VortexRe1000::primary_center_y};
    }
}

}
