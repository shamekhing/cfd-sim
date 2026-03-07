#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "../core/Field.hpp"

namespace io {

/**
 * @brief RGB Color representation
 */
struct RGB {
    float r, g, b;
    
    RGB() : r(0), g(0), b(0) {}
    RGB(float r_, float g_, float b_) : r(r_), g(g_), b(b_) {}
};

/**
 * @brief Color scheme enumeration
 */
enum class ColorScheme {
    BLUE_GREEN_RED,  // Blue (low) → Green (mid) → Red (high)
    VIRIDIS,         // Perceptually uniform: Blue → Green → Yellow
    HOT,             // Black → Red → Yellow → White
    COOL             // Cyan → Magenta
};

/**
 * @brief Color map utility for scalar field visualization
 * 
 * Generates RGB lookup tables and normalizes scalar fields to [0,1] range
 * for visualization in VTK format with custom color schemes.
 */
class ColorMap {
public:

    /**
     * @brief Map a normalized value [0,1] to RGB color
     * 
     * @param normalized_value Value in [0, 1] range
     * @param scheme Color scheme to use
     * @return RGB color triplet
     */
    static RGB mapValue(float normalized_value, ColorScheme scheme = ColorScheme::BLUE_GREEN_RED) {
        // Clamp to [0, 1]
        normalized_value = std::max(0.0f, std::min(1.0f, normalized_value));
        
        switch (scheme) {
            case ColorScheme::BLUE_GREEN_RED:
                return blueGreenRed(normalized_value);
            case ColorScheme::VIRIDIS:
                return viridis(normalized_value);
            case ColorScheme::HOT:
                return hot(normalized_value);
            case ColorScheme::COOL:
                return cool(normalized_value);
            default:
                return RGB(0.5f, 0.5f, 0.5f);  // Gray fallback
        }
    }

    /**
     * @brief Generate a lookup table with 256 colors
     * 
     * @param scheme Color scheme to use
     * @return Vector of 256 RGB colors
     */
    static std::vector<RGB> generateLookupTable(ColorScheme scheme = ColorScheme::BLUE_GREEN_RED) {
        std::vector<RGB> lut(256);
        for (int i = 0; i < 256; ++i) {
            float normalized = static_cast<float>(i) / 255.0f;
            lut[i] = mapValue(normalized, scheme);
        }
        return lut;
    }

    /**
     * @brief Find min and max values in a field
     * 
     * @param field Input scalar field
     * @return Pair of (min_value, max_value)
     */
    static std::pair<double, double> getFieldRange(const core::Field& field) {
        if (field.size() == 0) {
            return {0.0, 1.0};
        }
        
        double min_val = field[0];
        double max_val = field[0];
        
        for (int idx = 0; idx < field.size(); ++idx) {
            min_val = std::min(min_val, static_cast<double>(field[idx]));
            max_val = std::max(max_val, static_cast<double>(field[idx]));
        }
        
        // Avoid zero-range
        if (std::abs(max_val - min_val) < 1e-14) {
            min_val -= 0.5;
            max_val += 0.5;
        }
        
        return {min_val, max_val};
    }

    /**
     * @brief Normalize a field to [0, 1] range and map to RGB indices [0, 255]
     * 
     * @param field Input scalar field
     * @param min_val Minimum value for normalization
     * @param max_val Maximum value for normalization
     * @return Vector of indices into lookup table (0-255)
     */
    static std::vector<unsigned char> fieldToColorIndices(
        const core::Field& field,
        double min_val,
        double max_val)
    {
        std::vector<unsigned char> indices(field.size());
        
        double range = max_val - min_val;
        if (std::abs(range) < 1e-14) {
            std::fill(indices.begin(), indices.end(), 127);  // Gray
            return indices;
        }
        
        for (int idx = 0; idx < field.size(); ++idx) {
            double normalized = (static_cast<double>(field[idx]) - min_val) / range;
            normalized = std::max(0.0, std::min(1.0, normalized));
            
            unsigned char index = static_cast<unsigned char>(normalized * 255.0);
            indices[idx] = index;
        }
        
        return indices;
    }

private:

    /**
     * @brief Blue → Green → Red color mapping
     * 
     * v=0.0: (0, 0, 1) blue
     * v=0.5: (0, 1, 0) green
     * v=1.0: (1, 0, 0) red
     */
    static RGB blueGreenRed(float v) {
        if (v < 0.5f) {
            // Blue to Green: v ∈ [0, 0.5]
            float t = v * 2.0f;  // Rescale to [0, 1]
            return RGB(0.0f, t, 1.0f - t);
        } else {
            // Green to Red: v ∈ [0.5, 1.0]
            float t = (v - 0.5f) * 2.0f;  // Rescale to [0, 1]
            return RGB(t, 1.0f - t, 0.0f);
        }
    }

    /**
     * @brief Viridis-inspired color mapping (Blue → Green → Yellow)
     * 
     * Approximation of perceptually uniform colormap
     */
    static RGB viridis(float v) {
        // Viridis-like colormap (simplified)
        if (v < 0.25f) {
            float t = v / 0.25f;
            return RGB(0.267f + t * 0.093f, 0.004f + t * 0.329f, 0.329f + t * 0.340f);
        } else if (v < 0.5f) {
            float t = (v - 0.25f) / 0.25f;
            return RGB(0.360f + t * 0.310f, 0.333f + t * 0.338f, 0.669f - t * 0.380f);
        } else if (v < 0.75f) {
            float t = (v - 0.5f) / 0.25f;
            return RGB(0.670f + t * 0.265f, 0.671f + t * 0.250f, 0.289f - t * 0.289f);
        } else {
            float t = (v - 0.75f) / 0.25f;
            return RGB(0.935f + t * 0.065f, 0.921f + t * 0.079f, 0.0f);
        }
    }

    /**
     * @brief Hot color mapping (Black → Red → Yellow → White)
     * 
     * Good for highlighting hot regions
     */
    static RGB hot(float v) {
        if (v < 1.0f/3.0f) {
            float t = v / (1.0f/3.0f);
            return RGB(t, 0.0f, 0.0f);
        } else if (v < 2.0f/3.0f) {
            float t = (v - 1.0f/3.0f) / (1.0f/3.0f);
            return RGB(1.0f, t, 0.0f);
        } else {
            float t = (v - 2.0f/3.0f) / (1.0f/3.0f);
            return RGB(1.0f, 1.0f, t);
        }
    }

    /**
     * @brief Cool color mapping (Cyan → Magenta)
     * 
     * Complementary scheme for pressure/divergence
     */
    static RGB cool(float v) {
        return RGB(v, 0.0f, 1.0f - v);
    }
};

}
