#pragma once

#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "../include/assembly/GenericScalarAssembler.hpp"
#include "../config/ScalarConfig.hpp"
#include "../include/linear/BiCGSTAB.hpp"

namespace physics {

/**
 * Generic scalar transport equation solver
 * 
 * Solves: ∂φ/∂t + u·∇φ = ∇·(Γ∇φ) + S_φ
 * 
 * With optional:
 * - Viscous dissipation (for energy)
 * - Volumetric source terms
 * - Multiple convection schemes
 */
class ScalarTransport {
public:
	/**
	 * Solve scalar transport equation for one time step
	 *
	 * @param mesh          Structured mesh
	 * @param phi           Scalar field (temperature, concentration, etc.)
	 * @param u, v          Velocity components
	 * @param cfg           Scalar configuration (diffusivity, scheme, etc.)
	 * @param dt            Time step
	 * @param nu            Kinematic viscosity (for dissipation term)
	 * @param source_field  Optional volumetric source (nullptr = no source)
	 * @return              Number of linear solver iterations
	 */
	static int solveTransient(
		core::Field& phi,
		const core::Mesh& mesh,
		const core::Field& u,
		const core::Field& v,
		const config::ScalarConfig& cfg,
		linear::Scalar dt,
		linear::Scalar nu = 0.01,
		const core::Field* source_field = nullptr,
		int max_iter = 10);

	/**
	 * Apply Dirichlet boundary condition
	 */
	static void applyBC_Dirichlet(
		core::Field& phi,
		linear::Scalar bc_value = 0.0);
};

}  // namespace physics