#pragma once

#include "../core/Mesh.hpp"
#include "../core/Field.hpp"
#include "../core/Indexing.hpp"
#include "../linear/Vector.hpp"
#include "../linear/Matrix.hpp"
#include "../config/ScalarConfig.hpp"

namespace assembly {

// Forward declaration for convection scheme
enum class ConvectionScheme {
    Central,   // 2nd-order, non-monotonic
    Upwind,    // 1st-order, monotonic
    QUICK,     // 3rd-order, some oscillation
    MUSCL,     // 2nd-order TVD
    VanLeer    // 2nd-order TVD, smooth limiter
};

struct ScalarSystem {
    linear::SparseMatrix A;
    linear::Vector b;
};

class GenericScalarAssembler {
public:
    /**
     * Assemble generic scalar transport equation:
     * ∂φ/∂t + u·∇φ = ∇·(Γ∇φ) + S_φ
     *
     * For energy (is_energy=true):
     *   Γ = α (thermal diffusivity)
     *   S_φ includes viscous dissipation: Φ = 2μ(S·S) where S = (∇u + ∇u^T)/2
     *
     * @param mesh             Structured Cartesian mesh
     * @param phi              Scalar field (φ or T)
     * @param u, v             Velocity components
     * @param gamma            Diffusivity (Γ)
     * @param dt               Time step (for unsteady terms)
     * @param phi_old          Previous time step scalar (for ∂φ/∂t term)
     * @param scheme           Convection discretization scheme
     * @param source_field     Optional volumetric source S_φ (nullptr = zero source)
     * @param is_energy        If true, include viscous dissipation in source
     * @return                 Matrix A and RHS vector b
     */
    static ScalarSystem assemble(
        const core::Mesh& mesh,
        const core::Field& phi,
        const core::Field& u,
        const core::Field& v,
        double gamma,
        double dt,
        const core::Field& phi_old,
        ConvectionScheme scheme = ConvectionScheme::Upwind,
        const core::Field* source_field = nullptr,
        bool is_energy = false,
        double nu = 0.01);  // kinematic viscosity (for dissipation)

    /**
     * Compute viscous dissipation term Φ = 2μS·S
     * where S_ij = (∂u_i/∂x_j + ∂u_j/∂x_i)/2
     *
     * For 2D:
     * Φ = 2μ[ (∂u/∂x)² + (∂v/∂y)² + ((∂u/∂y + ∂v/∂x)/2)² ]
     *
     * @param mesh    Mesh
     * @param u, v    Velocity components
     * @param nu      Kinematic viscosity
     * @return        Dissipation field at each cell
     */
    static core::Field computeViscousDissipation(
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v,
        double nu);

private:
    // Helper: interpolate scalar to face using specified scheme
    static double interpolateToFace(
        double phi_far_upstream,
        double phi_upstream,    // u-side for u-face, v-side for v-face
        double phi_center,
        double phi_downstream,  // e-side for u-face, n-side for v-face
        double velocity,        // u or v component
        double grid_spacing,    // dx or dy
        ConvectionScheme scheme);

    // Helper: compute upwind-biased interpolation
    static double upwindInterpolation(
        double phi_upstream,
        double phi_center,
        double velocity,
        double grid_spacing);

    // Helper: compute central difference interpolation
    static double centralInterpolation(
        double phi_upstream,
        double phi_center,
        double phi_downstream,
        double grid_spacing);

    // Helper: compute QUICK interpolation (3-point, requires far-upstream)
    static double quickInterpolation(
        double phi_far_upstream,
        double phi_upstream,
        double phi_center,
        double phi_downstream,
        double velocity,
        double grid_spacing);

    // Helper: compute van Leer limiter (for TVD schemes)
    static double vanLeerLimiter(double r);

    // Helper: MUSCL slope limiter
    static double musclLimiter(double r);

    // Helper: assemble with specified scheme
    static ScalarSystem assembleWithScheme(
        const core::Mesh& mesh,
        const core::Field& phi,
        const core::Field& u,
        const core::Field& v,
        double gamma,
        double dt,
        const core::Field& phi_old,
        ConvectionScheme scheme,
        const core::Field* source_field,
        bool is_energy,
        double nu);
};

}  // namespace assembly
