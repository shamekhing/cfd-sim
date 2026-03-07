#include "../include/physics/ScalarTransport.hpp"

namespace physics {

int ScalarTransport::solveTransient(
    core::Field& phi,
    const core::Mesh& mesh,
    const core::Field& u,
    const core::Field& v,
    const config::ScalarConfig& cfg,
    linear::Scalar dt,
    linear::Scalar nu,
    const core::Field* source_field,
    int max_iter)
{
    core::Field phi_old = phi;

    int total_iters = 0;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Assemble scalar equation
        auto system = assembly::GenericScalarAssembler::assemble(
            mesh,
            phi,
            u,
            v,
            cfg.diffusivity,
            dt,
            phi_old,
            static_cast<assembly::ConvectionScheme>(cfg.convection_scheme),
            source_field,
            cfg.is_energy,
            nu);

        // Solve using BiCGSTAB
        linear::BiCGSTABSolver solver;
        linear::Vector solution = solver.solve(system.A, system.b);
        
        // Update field
        phi.raw() = solution;

        // Stop if converged (optional relaxation could be added)
        total_iters += cfg.max_scalar_iterations;
    }

    return total_iters;
}

void ScalarTransport::applyBC_Dirichlet(
    core::Field& phi,
    linear::Scalar bc_value)
{
    (void)phi;
    (void)bc_value;
}

}  // namespace physics
