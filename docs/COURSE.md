# CFD Course: Theory, Learning Path & Foundations

**Course:** Numerical Methods in Fluid Dynamics (LVA # 360012) 
**Instructor:** Stefan Pirker, Department of Particulate Flow Modelling, JKU Linz 
**Implementation:** CFD Engine - Educational CUDA-accelerated Solver 

This document provides the complete educational foundation for learning CFD through this project, combining theory, guided learning progression, and code references.

**Navigation:**
- [Course Overview](#course-overview) — What is CFD, learning outcomes
- [Learning Path](#learning-path) — 7-level progression (40-50 hours)
- [Reference Frames](#reference-frames) — Lagrangian vs Eulerian (key concept)
- [Governing Equations](#governing-equations-summary) — Conservation laws
- [Numerical Methods](#numerical-methods-foundations) — Discretization principles
- [Turbulence Introduction](#turbulence-modeling-introduction) — k-ε basics

**For other documentation:**
- Code implementation → See [TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md)
- Building & running → See [USER_GUIDE.md](docs/USER_GUIDE.md)
- Project overview → See [README.md](README.md)

---

## Course Overview

### What is CFD?

Computational Fluid Dynamics solves the equations governing fluid flow using numerical methods on computers. It bridges theory (governing equations) and experiments by:

1. **Discretizing** continuous equations into algebraic systems
2. **Solving** large linear/nonlinear systems efficiently
3. **Visualizing** and validating results against benchmarks

### Course Structure

This project implements concepts from the course in progressive difficulty:

```
Foundation
 ↓
Governing Equations (Mass, Momentum, Energy, Species conservation)
 ↓
Numerical Methods (Finite Volume, Discretization schemes)
 ↓
Physics Modules (Staggered grids, SIMPLE algorithm, Transport)
 ↓
Real Applications (Lid-driven cavity, Energy transport, Turbulence)
```

### Learning Outcomes

After working through this project, you will:

- [Complete] Understand **finite volume discretization** principles
- [Complete] Implement **SIMPLE pressure-velocity coupling** algorithm
- [Complete] Solve **Navier-Stokes equations** on staggered grids
- [Complete] Validate solutions against **benchmark problems**
- [Complete] Extend solvers with **new physics modules** (turbulence, multiphase)
- [Complete] Accelerate computations using **CUDA**
- [Complete] Generate animations with **parallel I/O** (async file writing)

### Project Status

[Complete] **Phase 1 Complete:**
- Staggered grid SIMPLE solver
- Lid-driven cavity benchmark validation (Re=100)
- Energy transport with MMS validation
- VTK/VTS/PVD animation export
- Async parallel I/O (20-50% speedup)
- Animation resume with field state loading

[Complete] **Phase 2 Complete (v2.0.0):**
- k-ε turbulence model with wall functions
- Wall y+ diagnostics for mesh validation
- High-Reynolds test cases (cavity, channel)
- Source term support in SIMPLE solver

Plan **Roadmap (Phase 3+):**
- Advanced turbulence models (k-ω SST, Spalart-Allmaras)
- Large Eddy Simulation (LES)
- Multiphase flows (free surface, particulate)
- Advanced preconditioners
- GPU memory optimization

### Key Concepts at a Glance

| Concept | Slide | Code | Section |
|---------|-------|------|---------|
| **Reynolds Number** | Intro_marked | `config/SolverConfig` | [Overview](#course-overview) |
| **Reference Frames** | GovEqI_marked | `core/Mesh` | [Reference Frames](#reference-frames) |
| **Mass Conservation** | GovEqI_marked | `solver/SIMPLE` | [Governing Equations](#governing-equations-summary) |
| **Energy Equation** | lecture-3-marked | `physics/ScalarTransport` | [Governing Equations](#governing-equations-summary) |
| **Finite Volume** | FVII_unmarked | `numerics/` | [Numerical Methods](#numerical-methods-foundations) |
| **Species Transport** | SpecHeat_marked | `cases/EnergyTransport` | [Governing Equations](#governing-equations-summary) |
| **Turbulence k-ε** | TurbI/II/III | `physics/TurbulenceKEpsilon` | [Turbulence](#turbulence-modeling-introduction) |
| **SIMPLE Algorithm** | TurbI | `solver/StaggeredSIMPLE` | [Numerical Methods](#numerical-methods-foundations) |

### How to Use This Document

**Follow the learning path:**
```
Slide (theory)
 ↓
Course section (concepts + equations)
 ↓
Code examples (docs/TECHNICAL_REFERENCE.md)
 ↓
Validation (docs/USER_GUIDE.md)
```

**Each section includes:**
- Physical equations with LaTeX
- Conceptual explanations
- References to implementation (see docs/TECHNICAL_REFERENCE.md)
- Links to validation cases (see docs/USER_GUIDE.md)
- Pointers to lecture slides

---

## Learning Path

This learning path guides you through the course materials → documentation → working code, following the lecture progression.

### Time Estimates

- **Level 1-3 (Foundations):** 6-9 hours
- **Level 4-5 (Discretization):** 8-10 hours 
- **Level 6 (SIMPLE):** 5-6 hours
- **Level 7 (Turbulence):** 8-10 hours
- **Total Core Path:** ~40-50 hours

---

### Level 1: Introduction & Foundations

**Lecture Lecture:** `Intro_marked.pdf` (31 slides) 
**⏱ Time:** 1-2 hours

**Learning objectives:**
- Understand scope of computational fluid dynamics
- Learn Reynolds number and flow regimes
- Meet the project structure

**Your journey:**

```
Read: Slides Intro_marked.pdf
 ↓
Read: Course Overview (above)
 ↓
Run: cd build && ./cfd --help
 ↓
See: All available flags (--animate, --name, --numTimeSteps)
```

**Key takeaways:**
- CFD bridges theory/experiments via numerical discretization
- Reynolds number determines flow character (laminar/turbulent)
- This project implements SIMPLE algorithm for Navier-Stokes
- CUDA acceleration available for larger problems

---

### Level 2: Governing Equations I

**Lecture Lecture:** `GovEqI_marked.pdf` (21 slides) 
** Documentation:** [Governing Equations](#governing-equations-summary), [Reference Frames](#reference-frames) 
**⏱ Time:** 3-4 hours

**Topics covered:**
- Mathematical conventions (tensors, vectors, scalars)
- **Lagrangian vs Eulerian reference frames** ← KEY CONCEPT
- Mass conservation (differential & integral forms)

**Your journey:**

```
Read: Slides GovEqI_marked.pdf
 ↓
Read: Reference Frames section (below)
 ↓
Read: Governing Equations Summary section (below)
 ↓
Study: Code at core/Mesh.hpp (grid discretization)
 solver/SIMPLE.hpp (mass equation assembly)
 ↓
Run: ./cfd --animate --numTimeSteps 10
 (Watch first 10 frames of lid cavity)
```

**Key concepts:**
- **Eulerian frame** = fixed grid (what we implement)
- **Lagrangian frame** = moving with fluid particles (theory)
- Conservation of mass: $\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) = 0$
- Staggered grid stores velocities at face centers, pressure at cell centers

**Code pointers:**
- `core/Mesh.hpp` — Grid discretization
- `core/Field.hpp` — Cell-centered fields
- `core/UField.hpp`, `core/VField.hpp` — Face-centered velocities
- `solver/StaggeredSIMPLE.hpp` — Mass & momentum equations

---

### Level 3: Governing Equations II (Energy & Species)

**Lecture Lecture:** `lecture-3-marked.pdf` (Energy confusion resolution) 
** Documentation:** [Governing Equations](#governing-equations-summary) 
**⏱ Time:** 2-3 hours

**Topics covered:**
- Energy conservation equation
- Species conservation (generic transport)
- Clarification of energy vs enthalpy

**Your journey:**

```
Read: Slides lecture-3-marked.pdf (focus on "Energy conservation")
 ↓
Read: Governing Equations Summary (below)
 ↓
Study: Code at physics/ScalarTransport.hpp
 cases/EnergyTransport.hpp
 ↓
Run: ./cfd --case energy_transport (run energy MMS test)
```

**Key equations:**
- Energy: $\frac{\partial (\rho T)}{\partial t} + \nabla \cdot (\rho \vec{u} T) = \nabla \cdot (\alpha \nabla T) + Q$
- Generic scalar: $\frac{\partial (\rho \phi)}{\partial t} + \nabla \cdot (\rho \vec{u} \phi) = \nabla \cdot (\Gamma \nabla \phi) + S_\phi$

**Code pointers:**
- `physics/ScalarTransport.hpp` — Scalar conservation
- `cases/EnergyTransport.hpp` — Energy transport case
- `numerics/ConvectionSchemes.hpp` — Convection discretization
- `numerics/Laplacian.hpp` — Diffusion discretization

---

### Level 4: Finite Volume Discretization

**Lecture Lecture:** `FVII_unmarked.pdf` (Finite Volume Part 7) 
** Documentation:** [Numerical Methods](#numerical-methods-foundations) 
**⏱ Time:** 4-5 hours

**Topics covered:**
- Finite volume method principles
- Convection-diffusion equation discretization
- Numerical vs physical convergence
- Upwinding and stability

**Your journey:**

```
Read: Slides FVII_unmarked.pdf
 ↓
Read: Numerical Methods Foundations (below)
 ↓
Study: Code at numerics/ConvectionSchemes.hpp
 numerics/Laplacian.hpp
 numerics/Gradient.hpp
 numerics/Divergence.hpp
 ↓
See: docs/TECHNICAL_REFERENCE.md for implementation details
```

**Key principles:**
- Discretize conservation laws over control volumes
- Approximate fluxes (convective, diffusive) at faces
- Choose discretization schemes (upwind, central, higher-order)
- Balance stability vs accuracy

**Code pointers:**
- `numerics/ConvectionSchemes.cpp` — Upwind, CDS, higher-order
- `numerics/Laplacian.hpp` — Discrete Laplacian (5-point stencil)
- `numerics/Gradient.hpp` — Gradient reconstruction
- `numerics/Divergence.hpp` — Divergence calculation
- `assembly/GenericScalarAssembler.hpp` — Linear system assembly

---

### Level 5: Species & Heat Transport

**Lecture Lecture:** `SpecHeat_marked.pdf` (33 slides) 
** Documentation:** [Governing Equations](#governing-equations-summary) 
**⏱ Time:** 3-4 hours

**Topics covered:**
- Species transport equations
- Heat transfer formulations
- Data-driven approaches

**Your journey:**

```
Read: Slides SpecHeat_marked.pdf
 ↓
Read: Governing Equations Summary (below)
 ↓
Study: Code at cases/EnergyTransport.hpp
 ↓
Run: ./cfd --case energy_transport
 Compare results vs analytical (MMS convergence)
```

**Key concepts:**
- Scalar transport extends momentum solver
- Heat/species differ only in transport coefficient
- Manufactured Solution Method (MMS) validates discretization
- Convergence: numerical error should decrease as $O(h^2)$

**Code pointers:**
- `physics/ScalarTransport.hpp` — Generic scalar solver
- `cases/EnergyTransport.hpp` — Energy transport test case

---

### Level 6: Pressure-Velocity Coupling (SIMPLE)

**Lecture Lecture:** `TurbI.pdf` (wrap-up review) 
** Documentation:** [Numerical Methods](#numerical-methods-foundations) 
**⏱ Time:** 5-6 hours

**Topics covered:**
- SIMPLE algorithm steps
- Pressure correction equation
- Checkerboard pressure oscillation prevention
- Staggered grids

**Your journey:**

```
Read: Slides TurbI.pdf (sections on Navier-Stokes coupling)
 ↓
Read: Numerical Methods Foundations (below)
 ↓
Study: docs/TECHNICAL_REFERENCE.md (SIMPLE implementation section)
 ↓
Study: Code at solver/StaggeredSIMPLE.hpp
 solver/SIMPLE.hpp
 ↓
Run: ./cfd --animate --numTimeSteps 50 (lid cavity animation)
```

**Algorithm: SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)**

```
1. Guess initial pressure field p*
2. Solve momentum equations with p* → u*, v*
3. Compute mass imbalance (residual)
4. Solve pressure correction equation → p'
5. Correct: p_new = p* + α_p * p'
6. Correct: u_new = u* + (∇p' correction)
7. Check convergence; if not, return to step 2
```

**Code pointers:**
- `solver/StaggeredSIMPLE.hpp` — Main SIMPLE solver
- `solver/SIMPLE.hpp` — Pure SIMPLE (steady-state)
- `assembly/MomentumAssembler.hpp` — Momentum equation LHS/RHS
- `assembly/PressureCorrectionAssembler.hpp` — Pressure correction

**For implementation details** → See [TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md#simple-algorithm)

---

### Level 7: Turbulence Modeling

**Lecture Lecture:** `TurbI.pdf`, `TurbII.pdf`, `TurbIII.pdf` (3 parts) 
** Documentation:** [Turbulence Introduction](#turbulence-modeling-introduction) 
**⏱ Time:** 8-10 hours

**Topics covered:**
- Reynolds decomposition
- Reynolds stress closure problem
- k-ε turbulence model
- Wall functions

**Your journey:**

```
Read: Slides TurbI/II/III.pdf
 ↓
Read: Turbulence Modeling Introduction (below)
 ↓
Study: docs/TECHNICAL_REFERENCE.md (turbulence implementation)
 ↓
Study: Code at physics/TurbulenceKEpsilon.hpp
 ↓
Run: ./cfd --case high_re_cavity --re 10000
 ./cfd --case channel_flow
 ↓
Validate: docs/USER_GUIDE.md (validation test cases)
```

**Key concepts:**
- Reynolds decomposition: $u = \bar{u} + u'$ (mean + fluctuation)
- Turbulent viscosity hypothesis: $\mu_t = C_\mu \rho \frac{k^2}{\epsilon}$
- k-ε model solves two transport equations (kinetic energy, dissipation)
- Wall functions bridge viscous sublayer without resolving it

**Code pointers:**
- `physics/TurbulenceKEpsilon.hpp` — k-ε model implementation
- Wall function implementation (iterative friction velocity solver)
- y+ diagnostics for mesh validation

**For complete implementation** → See [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md#turbulence-implementation) 
**For running turbulent cases** → See [USER_GUIDE.md](USER_GUIDE.md#turbulent-validation)

---

### Quick Reference: Slide ↔ Doc ↔ Code Map

| # | Lecture Slide | Topic | Course Section | Key Code |
|---|---|---|---|---|
| 1 | `Intro_marked.pdf` | CFD fundamentals | [Overview](#course-overview) | — |
| 2 | `GovEqI_marked.pdf` | Mass conservation | [Governing Equations](#governing-equations-summary) | core/, solver/ |
| 3 | `lecture-3-marked.pdf` | Energy equation | [Governing Equations](#governing-equations-summary) | physics/ScalarTransport |
| 4 | `FVII_unmarked.pdf` | Finite volume | [Numerical Methods](#numerical-methods-foundations) | numerics/, assembly/ |
| 5 | `SpecHeat_marked.pdf` | Species/heat | [Governing Equations](#governing-equations-summary) | cases/EnergyTransport |
| 6 | `TurbI.pdf` | SIMPLE coupling | [Numerical Methods](#numerical-methods-foundations) | solver/StaggeredSIMPLE |
| 7 | `TurbII/III.pdf` | Turbulence k-ε | [Turbulence](#turbulence-modeling-introduction) | physics/TurbulenceKEpsilon |

---

## Reference Frames

**Lecture:** `GovEqI_marked.pdf` (Key concept from slides)

### Overview

The **reference frame** is how we describe fluid motion. The choice fundamentally affects how we write and solve governing equations.

#### Two Perspectives

```
LAGRANGIAN FRAME EULERIAN FRAME
Follow fluid particles Fixed observation grid
(like a swimmer with the flow) (like a bridge watching a river)

Tracks particle paths Observes flow at fixed points
```

---

### Lagrangian Description

**Idea:** Label a fluid particle and follow it through space/time.

**Position:** $\vec{X}(a, t)$ where $a$ is initial particle label 
**Velocity:** $\vec{U} = \frac{D\vec{X}}{Dt}$ (following the particle) 
**Acceleration:** $\frac{D\vec{U}}{Dt}$ (rate of change along particle path)

#### Example
Imagine a smoke particle entering a jet:
- At time $t_0$: particle is at $(0, 0)$ with velocity $(1, 0)$
- At time $t_1$: particle is at $(1, 0.1)$ with velocity $(1.05, 0.02)$
- At time $t_2$: particle is at $(2.1, 0.15)$ with velocity $(1.1, 0.05)$

We track the **same particle** and observe its properties changing along its path.

#### Advantages
- [Complete] Natural for particle dynamics (sprays, sediment transport)
- [Complete] Conservation laws are simple (only time derivatives)
- [Complete] Lagrangian particle tracking used in multiphase flows

#### Disadvantages
- [Not Implemented] Expensive computationally (track millions of particles)
- [Not Implemented] Difficult to implement fixed boundary conditions
- [Not Implemented] Grid distortion problem (Lagrangian elements get tangled)

---

### Eulerian Description

**Idea:** Watch a fixed point in space as particles pass through it.

**Position:** $\vec{x}$ (fixed observation point) 
**Velocity field:** $\vec{u}(\vec{x}, t)$ (velocity at a point) 
**Acceleration:** $\frac{D\vec{u}}{Dt} = \frac{\partial \vec{u}}{\partial t} + (\vec{u} \cdot \nabla) \vec{u}$

#### Decomposition of Material Derivative

The **material derivative** (following a particle) splits into two parts:

$$\frac{D\phi}{Dt} = \frac{\partial \phi}{\partial t} + (\vec{u} \cdot \nabla) \phi$$

Where:
- $\frac{\partial \phi}{\partial t}$ = **Local time change** (at fixed point)
- $(\vec{u} \cdot \nabla) \phi$ = **Convection** (change due to particle movement)

#### Example

Temperature at a point near a heater:
- If heater is hotter → $\frac{\partial T}{\partial t} > 0$ (local increase)
- If wind blows hot air away → we feel cooler despite heater (convection term)
- Total change experienced by air particle = both effects combined

#### Advantages
- [Complete] **Natural for CFD** (fixed computational grid)
- [Complete] Easy to implement boundary conditions
- [Complete] No mesh distortion (Eulerian elements stay regular)
- [Complete] Efficient for engineering problems (pipe flow, aerodynamics)

#### Disadvantages
- [Not Implemented] Convection terms harder to discretize (numerical schemes crucial)
- [Not Implemented] Mass conservation requires divergence-free velocity field
- [Not Implemented] Tracking particles requires post-processing

---

### Comparison: Navier-Stokes Equations

#### Lagrangian Form (Following a particle)

$$\rho \frac{D\vec{u}}{Dt} = -\nabla p + \mu \nabla^2 \vec{u} + \vec{f}_{\text{body}}$$

**Interpretation:** The acceleration of a fluid particle equals forces acting on it. 
**Advantage:** Simple left-hand side (only material derivative) 
**Problem:** Right-hand side depends on Eulerian gradients (mixed description)

#### Eulerian Form (Fixed grid)

$$\rho \left( \frac{\partial \vec{u}}{\partial t} + (\vec{u} \cdot \nabla) \vec{u} \right) = -\nabla p + \mu \nabla^2 \vec{u} + \vec{f}_{\text{body}}$$

**Interpretation:** 
- Local momentum change + Momentum convection = Pressure + Viscous forces + Body forces

**Advantage:** All terms expressed at fixed points (compatible with fixed grid) 
**This is what we solve computationally.**

---

### Why Eulerian for This Project?

Our CFD solver uses **Eulerian description** because:

1. **Fixed grid efficiency** — Computations on same mesh every iteration
2. **Staggered grid compatibility** — Mesh doesn't move; structure is stable
3. **Boundary conditions** — Fixed walls, inlets, outlets are natural
4. **Standard practice** — Most CFD codes (ANSYS Fluent, OpenFOAM, etc.) use Eulerian
5. **Pressure handling** — Eulerian form separates pressure (Poisson equation)

#### The Price We Pay

Eulerian description introduces **convection terms** $(\vec{u} \cdot \nabla) \phi$, which are:
- **Nonlinear** in momentum equations (strong coupling)
- **Numerically tricky** (upwinding needed to prevent oscillations)
- **Expensive** (iteration required to converge)

See [Numerical Methods](#numerical-methods-foundations) for how we handle convection.

---

### Code Implementation

#### Grid (Eulerian)
```cpp
// Fixed mesh — Eulerian description
core::Mesh mesh(Nx, Ny, Lx, Ly); // Grid doesn't change during simulation
double x = mesh.x(i); // Position of cell i (fixed)
double u_at_x = u(i, j); // Velocity at that position
```

#### Mass Conservation (Eulerian form)
```cpp
// ∂ρ/∂t + ∇·(ρu) = 0

for (int i = 0; i < Nx; ++i) {
 for (int j = 0; j < Ny; ++j) {
 // Local change + Convection
 double local_change = du_dt(i, j);
 double convection = (u(i+1, j) - u(i-1, j)) / (2*dx);
 }
}
```

#### Material Derivative Implemented
```cpp
// D(u)/Dt = ∂u/∂t + (u·∇)u

auto material_accel = [&](int i, int j) {
 double local = u_new(i, j) - u_old(i, j); // ∂u/∂t
 double convection = u(i, j) * grad_u_x(i, j) 
 + v(i, j) * grad_u_y(i, j); // (u·∇)u
 return (local + dt * convection) / dt;
};
```

---

### Physical Insight

**When Convection Dominates (High Reynolds):**
- Particles carry momentum far from origin
- Eulerian equations need upwinding stabilization
- Material derivative term dominance

**When Diffusion Dominates (Low Reynolds):**
- Viscous forces spread momentum locally
- Eulerian discretization less problematic
- Laplacian term dominance

**In Our Lid-Driven Cavity (Re=100):**
- Moderate Reynolds number
- Both convection and diffusion important
- SIMPLE algorithm converges in ~5 iterations per timestep

---

### Summary

| Aspect | Lagrangian | Eulerian |
|--------|-----------|----------|
| **Follow** | Particle | Point in space |
| **Grid** | Moves with fluid | Fixed |
| **d/dt** | Material derivative only | Local + Convection |
| **Navier-Stokes LHS** | $\rho D\vec{u}/Dt$ | $\rho[\partial\vec{u}/\partial t + (\vec{u} \cdot \nabla) \vec{u}]$ |
| **Best for** | Particles, droplets | Pipe flow, aerodynamics |
| **Implementation** | Particle tracking | Fixed mesh CFD |
| **Used in this project** | [Not Implemented] Not used | [Complete] **Yes** |

---

## Governing Equations Summary

**The physics foundation — conservation laws and their mathematical forms**

### Overview

Fluid flow is governed by fundamental conservation laws:

1. **Mass Conservation** (Continuity) — No matter is created or destroyed
2. **Momentum Conservation** (Navier-Stokes) — Newton's 2nd law applied to fluids
3. **Energy Conservation** — Energy balance including heat transfer
4. **Species Transport** — Generic conservation equation for scalars

### Lecture References

| Topic | Slide |
|-------|-------|
| Reference frames | `GovEqI_marked.pdf` |
| Mass conservation | `GovEqI_marked.pdf` |
| Momentum equations | `GovEqI_marked.pdf` |
| Energy equation | `lecture-3-marked.pdf` |
| Species transport | `SpecHeat_marked.pdf` |

---

### Continuity (Mass Conservation)

**Differential form:**

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) = 0$$

**For incompressible flow** ($\rho =$ constant):

$$\nabla \cdot \vec{u} = 0$$

**Physical meaning:** Flow rate in = Flow rate out (no accumulation)

**Implementation:** `solver/StaggeredSIMPLE.hpp` enforces mass balance through pressure correction equation

---

### Momentum (Navier-Stokes)

**General form:**

$$\rho \frac{D\vec{u}}{Dt} = -\nabla p + \mu \nabla^2 \vec{u} + \vec{f}$$

**Expanded (Eulerian):**

$$\rho \left( \frac{\partial \vec{u}}{\partial t} + (\vec{u} \cdot \nabla) \vec{u} \right) = -\nabla p + \mu \nabla^2 \vec{u} + \vec{f}$$

**Physical meaning:**
- Left side: Momentum change (local + convection)
- Right side: Pressure forces + Viscous forces + Body forces

**Implementation:** `solver/StaggeredSIMPLE.hpp` solves momentum equations with SIMPLE coupling

---

### Energy (Temperature)

**General form:**

$$\rho c_p \frac{D T}{Dt} = \nabla \cdot (k \nabla T) + Q$$

**Expanded:**

$$\frac{\partial (\rho c_p T)}{\partial t} + \nabla \cdot (\rho c_p \vec{u} T) = \nabla \cdot (k \nabla T) + Q$$

**Physical meaning:**
- Left side: Energy storage + Convective transport
- Right side: Conductive heat transfer + Heat sources

**Implementation:** `physics/ScalarTransport.hpp` with thermal diffusivity $\alpha = k/(\rho c_p)$

---

### Species (Generic Scalar)

**General form:**

$$\frac{\partial (\rho \phi)}{\partial t} + \nabla \cdot (\rho \vec{u} \phi) = \nabla \cdot (\Gamma \nabla \phi) + S_\phi$$

Where:
- $\phi$ = scalar field (concentration, temperature, turbulent kinetic energy)
- $\Gamma$ = diffusivity coefficient
- $S_\phi$ = source term
- $\vec{u}$ = velocity field (from SIMPLE solver)

**Physical meaning:**
- Left side: Property accumulation + Convective transport
- Right side: Diffusive transport + Sources/sinks

**Implementation:** `physics/ScalarTransport.hpp` provides generic scalar solver framework

---

### Implementation Status

| Equation | Status | Code Location | Used In |
|----------|--------|---------------|---------|
| Mass conservation | [Complete] Complete | `solver/StaggeredSIMPLE.hpp` | All cases |
| Momentum (Navier-Stokes) | [Complete] Complete | `solver/StaggeredSIMPLE.hpp` | All cases |
| Energy transport | [Complete] Complete | `physics/ScalarTransport.hpp` | Energy transport case |
| Generic species | [Complete] Complete | `physics/ScalarTransport.hpp` | Custom scalars |
| Turbulence k-ε | [Complete] Complete | `physics/TurbulenceKEpsilon.hpp` | High-Re cases |

**For implementation details** → See [TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md)

---

## Numerical Methods Foundations

**How we discretize equations for computers**

### Overview

Converting continuous partial differential equations into discrete algebraic systems requires:

1. **Space discretization** — Divide domain into control volumes (finite volume method)
2. **Time discretization** — Approximate time derivatives (explicit/implicit schemes)
3. **Flux approximation** — Evaluate convective and diffusive fluxes at cell faces
4. **Linear system solution** — Solve resulting matrix equations efficiently

### Lecture Reference

| Topic | Slide |
|-------|-------|
| Finite Volume | `FVII_unmarked.pdf` |
| Convection-Diffusion | `FVII_unmarked.pdf` |
| Convergence & Stability | `FVII_unmarked.pdf` |

---

### Finite Volume Method

**Principle:** Integrate conservation laws over control volumes

$$\int_{CV} \frac{\partial \phi}{\partial t} dV + \int_{\partial CV} \vec{u} \phi \cdot \vec{n} dA = \int_{\partial CV} \Gamma \nabla \phi \cdot \vec{n} dA + \int_{CV} S_\phi dV$$

**Discretization:**
- Volume integral → value at cell center × cell volume
- Surface integral → sum of fluxes at face centers
- Results in algebraic equations for cell values

**Key features:**
- Inherently conservative (fluxes match at interfaces)
- Flexible for complex geometries
- Natural for incompressible flow (staggered grids)

**Implementation:** See [TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md#core-infrastructure) for mesh and field structures

---

### Convection-Diffusion Balance

**Generic transport equation:**

$$\frac{\partial \phi}{\partial t} + \nabla \cdot (\vec{u} \phi) = \nabla \cdot (\Gamma \nabla \phi)$$

**Péclet number** determines dominant mechanism:

$$Pe = \frac{\text{Convection}}{\text{Diffusion}} = \frac{|u| L}{\Gamma}$$

**Discretization challenges:**
- **High Péclet (convection-dominated):** Need upwinding for stability
- **Low Péclet (diffusion-dominated):** Central differencing accurate
- **Moderate Péclet:** Hybrid schemes balance accuracy and stability

---

### Discretization Schemes

**For convection term** $\nabla \cdot (\vec{u} \phi)$:

| Scheme | Accuracy | Stability | When to Use |
|--------|----------|-----------|-------------|
| **Upwind** | 1st order | Very stable | High-Re flows, debugging |
| **Central Difference** | 2nd order | Conditionally stable | Low-Re flows |
| **QUICK** | 3rd order | Conditionally stable | Smooth solutions |
| **MUSCL** | TVD scheme | Stable | Shock capturing |
| **Van Leer** | TVD scheme | Stable | Discontinuities |

**For diffusion term** $\nabla \cdot (\Gamma \nabla \phi)$:

- **Central differencing** (2nd order accurate)
- **5-point stencil** in 2D (7-point in 3D)

**Implementation:** See `numerics/ConvectionSchemes.cpp` and [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md#discretization-layer) for complete details

---

### Staggered Grid

**Problem:** Pressure-velocity coupling can produce checkerboard oscillations on collocated grids

**Solution:** Store variables at different locations

```
+-------+-------+
| | |
| p u p | p = pressure (cell center)
| | | u = x-velocity (vertical face)
+---v---+---v---+ v = y-velocity (horizontal face)
| | |
| p u p |
| | |
+-------+-------+
```

**Advantages:**
- Natural velocity divergence computation
- Pressure gradient at face matches velocity location
- Prevents checkerboard pressure oscillations

**Implementation:** 
- Pressure: `core/Field.hpp` (cell-centered)
- X-velocity: `core/UField.hpp` (vertical face-centered)
- Y-velocity: `core/VField.hpp` (horizontal face-centered)

**For staggered grid details** → See [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md#staggered-grid-implementation)

---

### SIMPLE Algorithm Concept

**Problem:** Momentum and continuity equations are coupled (pressure appears in momentum, velocity in continuity)

**SIMPLE** = **S**emi-**I**mplicit **M**ethod for **P**ressure-**L**inked **E**quations

**Algorithm steps:**

```
1. Guess initial pressure field p*
2. Solve momentum equations with p* → get u*, v* (intermediate velocities)
3. Compute mass imbalance: ∇·u* ≠ 0
4. Solve pressure correction equation → get p' to enforce ∇·u = 0
5. Correct pressure: p_new = p* + α_p * p'
6. Correct velocities: u_new = u* - (correction term involving ∇p')
7. Check convergence; if not converged, return to step 2
```

**Key parameters:**
- **Underrelaxation factors** $\alpha_p$ (pressure), $\alpha_u$ (velocity)
- Control convergence rate and stability
- Typical values: $\alpha_p \approx 0.3$, $\alpha_u \approx 0.7$

**Physical insight:**
- Momentum equations predict velocity (but may violate mass conservation)
- Pressure correction enforces mass conservation
- Iteration converges to solution satisfying both momentum and continuity

**Implementation:** See [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md#simple-algorithm) for complete algorithm with quasi-transient stabilization

---

### Time Integration

**Explicit schemes:**
- Easy to implement
- Restricted by CFL stability condition: $\Delta t \leq \frac{\Delta x}{|u|}$
- Not used for incompressible (stiff pressure coupling)

**Implicit schemes:**
- Unconditionally stable (can use larger timesteps)
- Require solving linear systems each timestep
- Used in this project for robustness

**Quasi-transient SIMPLE:** (Phase 1 innovation)
- Add pseudo-time term to momentum equations
- Improves convergence without affecting steady-state
- Allows larger underrelaxation factors

---

### Implementation Status

| Method | Status | Code Location |
|--------|--------|---------------|
| Finite Volume discretization | [Complete] Complete | `numerics/` + `assembly/` |
| Upwind convection scheme | [Complete] Complete | `numerics/ConvectionSchemes.cpp` |
| 5-point Laplacian | [Complete] Complete | `numerics/Laplacian.hpp` |
| Staggered grid layout | [Complete] Complete | `core/UField.hpp`, `core/VField.hpp` |
| SIMPLE coupling | [Complete] Complete | `solver/StaggeredSIMPLE.hpp` |
| Quasi-transient stabilization | [Complete] Complete | `solver/StaggeredSIMPLE.hpp` |

**For complete implementation** → See [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md)

---

## Turbulence Modeling Introduction

**Specialized equations for turbulent flows**

### Overview

Turbulence is characterized by:
- **Chaotic** flow with random fluctuations
- **Multi-scale** structures (large eddies to small dissipation)
- **Enhanced mixing** (momentum, heat, mass transfer)
- **High Reynolds numbers** (Re > 2300 for pipe flow)

**Challenge:** Direct simulation (DNS) requires resolving all scales — computationally prohibitive for engineering problems

**Solution:** Turbulence models provide closure for time-averaged equations (RANS - Reynolds-Averaged Navier-Stokes)

### Lecture References

| Topic | Slide |
|-------|-------|
| Turbulence Intro | `TurbI.pdf` |
| Turbulence Modeling | `TurbII.pdf` |
| Advanced Models | `TurbIII.pdf` |

---

### Reynolds Decomposition

**Idea:** Decompose flow variables into mean and fluctuating parts

$$u = \bar{u} + u', \quad v = \bar{v} + v', \quad p = \bar{p} + p'$$

Where:
- $\bar{u}$ = time-averaged velocity
- $u'$ = fluctuating component (with $\overline{u'} = 0$)

**Time averaging:**

$$\bar{u} = \frac{1}{T} \int_0^T u(t) dt$$

**Reynolds-averaged Navier-Stokes (RANS):**

$$\rho \left( \frac{\partial \bar{u}_i}{\partial t} + \bar{u}_j \frac{\partial \bar{u}_i}{\partial x_j} \right) = -\frac{\partial \bar{p}}{\partial x_i} + \frac{\partial}{\partial x_j} \left[ \mu \frac{\partial \bar{u}_i}{\partial x_j} - \rho \overline{u'_i u'_j} \right]$$

**New term:** $-\rho \overline{u'_i u'_j}$ = **Reynolds stress tensor** (turbulent momentum transport)

**Closure problem:** We have more unknowns ($\overline{u'_i u'_j}$) than equations!

---

### Turbulent Viscosity Hypothesis (Boussinesq)

**Assumption:** Reynolds stresses behave like viscous stresses

$$-\rho \overline{u'_i u'_j} = \mu_t \left( \frac{\partial \bar{u}_i}{\partial x_j} + \frac{\partial \bar{u}_j}{\partial x_i} \right) - \frac{2}{3} \rho k \delta_{ij}$$

Where:
- $\mu_t$ = **turbulent (eddy) viscosity** (not a fluid property!)
- $k = \frac{1}{2} \overline{u'_i u'_i}$ = turbulent kinetic energy
- $\delta_{ij}$ = Kronecker delta

**Now need:** Model for $\mu_t$ (turbulence closure model)

---

### k-ε Turbulence Model

**Two-equation model:** Solves transport equations for $k$ (turbulent kinetic energy) and $\epsilon$ (dissipation rate)

**Turbulent kinetic energy:**

$$\frac{\partial (\rho k)}{\partial t} + \nabla \cdot (\rho \vec{u} k) = \nabla \cdot \left[ \left( \mu + \frac{\mu_t}{\sigma_k} \right) \nabla k \right] + P_k - \rho \epsilon$$

**Dissipation rate:**

$$\frac{\partial (\rho \epsilon)}{\partial t} + \nabla \cdot (\rho \vec{u} \epsilon) = \nabla \cdot \left[ \left( \mu + \frac{\mu_t}{\sigma_\epsilon} \right) \nabla \epsilon \right] + C_{1\epsilon} \frac{\epsilon}{k} P_k - C_{2\epsilon} \rho \frac{\epsilon^2}{k}$$

**Turbulent viscosity:**

$$\mu_t = C_\mu \rho \frac{k^2}{\epsilon}$$

**Constants:**
- $C_\mu = 0.09$, $C_{1\epsilon} = 1.44$, $C_{2\epsilon} = 1.92$
- $\sigma_k = 1.0$, $\sigma_\epsilon = 1.3$

**Production term:** $P_k = \mu_t S^2$ (turbulence generation from mean flow shear)

---

### Wall Functions

**Problem:** Near walls, turbulent length scale $\sim y$ (distance from wall) → requires very fine mesh

**Solution:** Bridge near-wall region with empirical **wall functions**

**Logarithmic law of the wall:**

$$u^+ = \frac{1}{\kappa} \ln(y^+) + B$$

Where:
- $u^+ = u / u_\tau$ = dimensionless velocity
- $y^+ = \rho u_\tau y / \mu$ = dimensionless wall distance
- $u_\tau = \sqrt{\tau_w / \rho}$ = friction velocity
- $\kappa = 0.41$, $B = 5.2$ (empirical constants)

**Wall function application:**

1. Compute wall shear stress: $\tau_w = \rho u_\tau^2$
2. Enforce logarithmic profile at first grid point (instead of resolving viscous sublayer)
3. Set k and ε boundary conditions based on wall equilibrium assumptions

**Requirement:** First grid point at $y^+ \approx 30$ (log-law region)

**Implementation:** `physics/TurbulenceKEpsilon.hpp::applyWallFunctions()` with iterative friction velocity solver

---

### When to Use Turbulence Models

| Reynolds Number | Regime | Model Choice |
|----------------|--------|--------------|
| Re < 2300 | Laminar | No turbulence model (just Navier-Stokes) |
| 2300 < Re < 10,000 | Transitional | k-ε with care (or transition model) |
| Re > 10,000 | Fully turbulent | k-ε with wall functions (this project) |
| Re >> 100,000 | High-Re engineering | k-ω SST (better near-wall) |

**Example cases:**
- Lid-driven cavity Re=100 → **Laminar** (no turbulence model)
- Lid-driven cavity Re=10,000 → **k-ε with wall functions**
- Channel flow Re=10,000 → **k-ε with wall functions**

---

### Implementation Status (Phase 2 Complete!)

| Component | Status | Code Location |
|-----------|--------|---------------|
| k-ε transport equations | [Complete] Complete | `physics/TurbulenceKEpsilon.hpp` |
| Wall functions | [Complete] Complete | `physics/TurbulenceKEpsilon.hpp::applyWallFunctions()` |
| y+ diagnostics | [Complete] Complete | `solver/ResidualMonitor.hpp::wallYPlus()` |
| High-Re test cases | [Complete] Complete | Cavity, channel flow |
| Production term | [Complete] Complete | $P_k$ from strain rate tensor |
| Iterative u_tau solver | [Complete] Complete | Newton iteration for friction velocity |

**For complete implementation** → See [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md#turbulence-implementation) 
**For running turbulent cases** → See [USER_GUIDE.md](USER_GUIDE.md#turbulent-validation)

---

### Future Turbulence Models (Roadmap)

Plan **Planned (Phase 3+):**
- **k-ω SST model** — Better near-wall treatment, no wall functions needed
- **Spalart-Allmaras** — One-equation model for aerodynamics
- **Large Eddy Simulation (LES)** — Resolve large eddies, model small scales
- **Time-averaging infrastructure** — Compute turbulent statistics from simulations

**Estimated effort:** 10-25 hours per model

---

## References

- **Course:** LVA # 360012 "Numerical Methods in Fluid Dynamics"
- **Instructor:** Stefan Pirker (stefan.pirker@jku.at)
- **Slides:** See `course/slides/` folder for lecture PDFs
- **Implementation:** C++17 with CUDA backend (optional)
- **Validation Data:** Ghia et al. lid-driven cavity, analytical MMS solutions

---

## See Also

- **[TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md)** — Complete code architecture and implementation details
- **[USER_GUIDE.md](USER_GUIDE.md)** — Build instructions, running simulations, validation
- **[README.md](README.md)** — Project overview and quick start
