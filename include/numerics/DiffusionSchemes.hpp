#pragma once

#include "../linear/Vector.hpp"
#include "../../include/core/Mesh.hpp"
#include "../../include/core/Field.hpp"

namespace numerics {

/**
 * @brief Diffusion discretization schemes
 * 
 * Currently empty - supports future expansion to:
 * - Upwind diffusion (1st order, stable)
 * - Central difference diffusion (2nd order, standard)
 * - Higher-order schemes (compact, ENO)
 */
class DiffusionSchemes {
public:
    // Future implementations
};

}  // namespace numerics
