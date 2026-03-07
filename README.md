# CFD Engine — Vortex Shedding around a Circular Cylinder

A modular Computational Fluid Dynamics solver in C++17 implementing the transient SIMPLE algorithm for incompressible Navier-Stokes equations on staggered MAC grids, with optional GPU acceleration via CUDA.

**Course:** LVA # 360012 "Numerical Methods in Fluid Dynamics" (JKU Linz) 
**Implementation:** Modern C++17, Eigen3 (CPU) / CUDA (GPU)

---

## Quick Start

```bash
cd build
cmake .. && make -j4
./cfd --help
./cfd --animate --numTimeSteps 50 # 50-frame animation
./cfd --case energy_transport # Validation test
```

**Full build instructions** → [USER_GUIDE.md](docs/USER_GUIDE.md#installation)

---

## Documentation Overview

**4 comprehensive markdown files (replaces 20+ previous files):**

### 1. [COURSE.md](COURSE.md) — Learning & Theory
- Course overview and learning outcomes
- 7-level learning path (40-50 hours)
- Reference frames (Lagrangian vs Eulerian)
- Governing equations (mass, momentum, energy, species)
- Numerical methods foundations (FVM, discretization, SIMPLE)
- Turbulence modeling introduction (k-ε, wall functions)

**Start here if:** You're learning CFD from scratch

### 2. [TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md) — Code Architecture
- 4-layer architecture design
- Core infrastructure (mesh, fields, staggered grids)
- Discretization operators and schemes
- SIMPLE algorithm details
- Turbulence implementation (k-ε, wall functions, y+ diagnostics)
- Assembly layers and solver innovations

**Start here if:** You're extending or debugging the code

### 3. ▶ [USER_GUIDE.md](docs/USER_GUIDE.md) — Building & Running
- Installation and build options
- CLI commands and configuration
- Output file formats (VTK, CSV)
- ParaView visualization setup
- Validation test cases (cavity, energy transport, turbulent flows)
- Troubleshooting common issues

**Start here if:** You want to run simulations

### 4. Plan This README — Project Overview

---

## Current Status

**Version:** 2.0.0 — **Phase 2 Complete** 
**Status:** Turbulent incompressible Navier-Stokes with k-ε RANS + wall functions

### Implemented Components

**Core**:
- Structured Cartesian mesh, scalar/vector fields, indexing
- Sparse matrix assembly via triplets
- Linear solvers: CG, BiCGSTAB, Direct (SparseLU)

**Incompressible Navier-Stokes**:
- SIMPLE pressure-velocity coupling (steady-state + transient modes)
- Staggered grid (MAC) discretization
- Bounded upwind convection (diagonal dominance)
- Quasi-transient stabilization for steady-state solves

**Scalar Transport** (Phase 1):
- Generic convection-diffusion equation
- Five convection schemes: Central, Upwind, QUICK, MUSCL, VanLeer
- Flux limiters (van Leer, MUSCL, Superbee)
- Scalar dissipation / source terms

**Turbulence Modeling** (NEW - Phase 2):
- Standard k-ε RANS model (Launder & Spalding)
- Logarithmic wall functions for high-Re flows (y+ > 11.25)
- Automatic wall treatment selection (Re-based)
- Wall y+ diagnostics and mesh validation
- Source term support in SIMPLE (channel flow, body forces)

**Validation**:
- Lid-driven cavity benchmark (Re=10, 100, 1000)
- High-Reynolds turbulent cavity (Re=1000+, k-ε model)
- Turbulent channel flow (Re_τ=5600, pressure gradient driven)
- Ghia et al. (1982) centerline L2 error comparisons
- Energy transport manufactured solution test (PASS: L2 error 0.00597)
- Wall y+ validation (log-layer verification)
- Mass residual L2 monitoring
- VTK and CSV output

### Running Cases

```
./cfd # Lid cavity Re=100 (default)
./cfd --re 10 # Lid cavity Re=10
./cfd --re 1000 # Lid cavity Re=1000
./cfd --case energy_transport # Energy transport test
./cfd --case high_re_cavity --re 1000 # Turbulent cavity with k-ε
./cfd --case channel_flow --re 5600 # Turbulent channel flow
./cfd --help # Show available options
```

### Phase 1 Test Results

**Lid-Driven Cavity (Re=100)**:
- SIMPLE iterations: 300 (converged)
- Mass residual L2: 1.94e-19
- Max velocity (interior): 0.796
- Ghia centerline L2 error (u): 0.344
- Ghia centerline L2 error (v): 0.226

**Energy Transport** (Manufactured Solution):
- Grid: 64×64, α=0.01, dt=1e-4, 100 time steps
- L2 error: 0.00597
- Relative error: 0.60%
- Status: **PASS**

### Next Milestone (Phase 3)

- Time-averaging infrastructure for turbulent statistics
- Additional turbulence models (k-ω SST, Spalart-Allmaras)
- Large Eddy Simulation (LES) framework
- Advanced time integration schemes (RK2, RK3-TVD, BDF2)
- Unsteady cylinder flow with vortex shedding

------------------------------------------------------------------------

## Flagship Case: Cylinder Vortex Shedding

**Source:** [`cases/main/CylinderMain.hpp`](cases/main/CylinderMain.hpp) 
**Config:** [`include/config/CaseConfigs.hpp`](include/config/CaseConfigs.hpp) → `CylinderFlowConfig`
**Sample Video:** [`cases/main/cyl_200x100re1000stp1500.avi`](cases/main/cyl_200x100re1000stp1500.avi)

### Problem Description

Flow of a viscous incompressible fluid past a circular cylinder placed at the centre of a rectangular domain. At Reynolds number Re = 100 the wake becomes unstable and sheds a periodic train of alternating vortices — the **Kármán vortex street**. This is the canonical unsteady laminar benchmark in computational fluid dynamics.

```
 U_inf → ← zero-gradient outlet
 ─────────────────────────────────────────────────────
 U_inf → ╭───╮
 U_inf → │ ● │ ← cylinder, centre (Lx/2, Ly/2)
 U_inf → ╰───╯ ↔ vortex shedding wake
 U_inf →
 ─────────────────────────────────────────────────────
 inlet outlet

 Domain : [0, Lx] × [0, Ly] default 40 × 24
 Cylinder diameter D = 4.0, radius r = 2.0
```

### Domain & Parameters

| Parameter | Symbol | Value | Notes |
|-----------|--------|-------|-------|
| Domain length | $L_x$ | 40.0 | 10 diameters |
| Domain height | $L_y$ | 24.0 | 6 diameters |
| Cylinder diameter | $D$ | 4.0 | reference length |
| Cylinder centre | $(c_x, c_y)$ | $(L_x/2,\; L_y/2)$ | centred cylinder |
| Free-stream velocity | $U_\infty$ | 1.0 | |
| Reynolds number | $Re$ | 100 | $Re = U_\infty D / \nu$ |
| Kinematic viscosity | $\nu$ | 0.04 | derived from Re |
| Default grid | $(N_x, N_y)$ | (120, 72) | square cells ($\Delta x=\Delta y=1/3$), $D/\Delta x\approx12$ |
| Time step | $\Delta t$ | 0.05 | CFL ≈ 0.15 on default grid |
| Default steps | — | 5 000 | covers ~42 shedding periods |
| Inner SIMPLE iters | — | 5 | robust transient SIMPLE convergence |

### Boundary Conditions

| Boundary | $u$ | $v$ | $p'$ |
|----------|-----|-----|------|
| Left (inlet) | $U_\infty\,r(t)$ (uniform ramp) | 0 | zero-gradient |
| Right (outlet) | $\partial u/\partial x = 0$ | 0 | 0 |
| Top / Bottom | $U_\infty$ (free-slip) | 0 | zero-gradient |
| Cylinder surface | 0 (IBM mask) | 0 (IBM mask) | 0 |

Startup ramp (first 2 s):

$$r(t)=\begin{cases}
	frac{1}{2}\left(1-\cos\left(\pi t/2\right)\right), & t<2 \\
1, & t\ge 2
\end{cases}$$

The case also applies a one-shot antisymmetric $v$ impulse at $t=0$ in the near wake to seed the unstable mode; sustained oscillation must persist long after this transient decays.

### Physics

**Governing equations** (incompressible Navier-Stokes):

$$\nabla \cdot \vec{u} = 0$$

$$\frac{\partial \vec{u}}{\partial t} + \vec{u} \cdot \nabla \vec{u} = -\nabla p + \nu \nabla^2 \vec{u}$$

**Expected flow regimes:**

| Re | Wake behaviour |
|----|---------------|
| < 5 | Fully attached, creeping flow |
| 5 – 40 | Steady symmetric recirculation pair |
| 40 – 47 | Steady asymmetric wake |
| **47 – 180** | **Periodic Kármán vortex shedding ← Re=100 is here** |
| 180 – 300 | Three-dimensional transition |
| > 300 | Turbulent wake |

### Numerical Method

The transient SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm on a **staggered MAC (Marker and Cell) grid**:

```
For each time step:
 u_old = u, v_old = v

 For it = 1 … itersPerTimeStep (= 5):
 1. Assemble u-momentum (implicit, with u_old/dt inertia)
 Override boundary rows → correct inlet / outlet / slip BCs
 Solve A_u · u = b_u

 2. Assemble v-momentum (implicit, with v_old/dt inertia)
 Solve A_v · v = b_v

 3. Apply IBM mask: zero all faces inside the cylinder

 4. Solve pressure-correction Poisson equation
 A_p · p' = divergence(u, v)

 5. Correct velocities:
 u(i,j) += (dy / a_e) · (p'_{i-1,j} - p'_{i,j})
 v(i,j) += (dx / a_n) · (p'_{i,j-1} - p'_{i,j})

 6. Update pressure: p += α_p · p'

 7. Re-enforce BCs + IBM mask

 Write animation frame (PVD/VTS)
 Compute Cd, Cl → forces.csv
 Estimate Strouhal → strouhal.csv
```

**Key choice: 2–5 inner SIMPLE iterations per time step (default 5).** 
Over-converging (e.g. 200 iterations at each step) drives the solution toward pseudo-steady behavior and can suppress physical unsteadiness.

### Symmetry-Breaking Strategy

The geometry is centred and symmetric. The unstable mode is seeded by a one-shot antisymmetric near-wake $v$ impulse at $t=0$, while the inlet itself remains uniform (with startup ramp). This avoids permanently biasing the wake direction.

### Immersed Boundary Method (IBM)

The cylinder is represented on the Cartesian grid via the **volume-penalisation / direct forcing** approach:

- After each momentum solve, any u-face whose staggered position $(i \Delta x,\, (j+\tfrac{1}{2})\Delta y)$ satisfies $d \leq r$ is set to zero.
- Similarly for v-faces at $((i+\tfrac{1}{2})\Delta x,\, j \Delta y)$.
- No body-fitted mesh required — the cylinder boundary is implicitly defined by which faces are zeroed.

Accuracy scales as $O(\Delta x)$; true Cd converges to the literature value as the mesh is refined (D/$\Delta x$ > 20 recommended for quantitative results).

### How to verify shedding is present

Use all three checks together:

1. **Flow-field signature (qualitative):** alternating vortices in the wake (staggered upper/lower structures), not a steady symmetric recirculation bubble.
2. **Lift signal (time series):** `Cl` must oscillate with repeated sign changes over many periods (not just a startup bump).
3. **Frequency metric (quantitative):** non-zero, stable Strouhal estimate from late-time `Cl` history.

Important visualization note: do **not** interpret `velocity_rgb Magnitude` as physical velocity; it is the magnitude of an RGB color vector. For physics checks, use scalar `velocity` magnitude (or vorticity).

### Validation Targets

| Quantity | This solver (default 120×72) | Literature | Reference |
|----------|-------------------|------------|-----------|
| Strouhal number $St = fD/U_\infty$ | typically 0.12–0.14 (late-time, current setup) | **0.164** | Williamson (1989) |
| Mean drag coefficient $\overline{C_d}$ | converges with mesh | **1.3 ± 0.05** | Williamson (1989) |
| Lift coefficient oscillation | periodic at $f_{shed}$ | periodic | — |
| Onset of shedding | after startup transient | — | — |

The Strouhal number is less grid-sensitive than drag, but accurate $\overline{C_d}$ and benchmark-level $St$ still require finer IBM resolution (target $D/\Delta x \gtrsim 20$).

### Force & Strouhal Computation

**Drag and lift coefficients** are computed by integrating pressure and approximate viscous stress over the one-cell-wide annular band just outside the cylinder surface:

$$C_d = \frac{F_x}{\frac{1}{2} \rho U_\infty^2 D}, \qquad C_l = \frac{F_y}{\frac{1}{2} \rho U_\infty^2 D}$$

Each surface cell is assigned an equal arc-length share $\Delta s = \pi D / N_{surface}$ to keep the total integral equal to the cylinder half-circumference regardless of grid resolution.

**Strouhal number** is detected online from the lift-coefficient history via zero-crossing period:

$$St = \frac{D}{U_\infty \cdot T_{shed}}, \qquad T_{shed} = \Delta t \cdot (k_n - k_{n-1})$$

where $k_n$ are successive Cl sign-change indices.

### Running the Case

```bash
# Quick smoke test (3 min on CPU)
./build_cuda/cfd --case cylinder_main --nx 64 --ny 64 --steps 500 --re 100

# Full validation run (recommended)
./build_cuda/cfd --case cylinder_main --nx 120 --ny 72 --re 100 --name flagship

# High-resolution quantitative Cd (slow)
./build_cuda/cfd --case cylinder_main --nx 256 --ny 128 --re 100 --name hires
```

**Resume support:** restart an interrupted run without losing frames:
```bash
# First run (interrupted at step 2000)
./build_cuda/cfd --case cylinder_main --nx 120 --ny 72 --name flagship

# Resume automatically picks up from step 2000
./build_cuda/cfd --case cylinder_main --nx 120 --ny 72 --name flagship
```

### Output Files

| File | Contents |
|------|----------|
| `output/flagship/vts/flagship_NNNNNNN.vts` | Per-frame velocity + pressure fields |
| `output/flagship/flagship.pvd` | ParaView time-series collection |
| `output/flagship/csv/forces.csv` | `step, time, Cd, Cl` per time step |
| `output/flagship/csv/strouhal.csv` | `step, time, f_shed, St` (when detected) |

### Results Media


**Cylinder snapshots by Reynolds number** (filename code uses scientific notation: `1e02` = Re=100):

| Reynolds number | Snapshot |
|---|---|
| Re = 10    | [`cases/main/results/flagship1e01.png`](cases/main/results/flagship1e01.png) |
| Re = 30    | [`cases/main/results/flagship3e01.png`](cases/main/results/flagship3e01.png) |
| Re = 60    | [`cases/main/results/flagship6e01.png`](cases/main/results/flagship6e01.png) |
| Re = 70    | [`cases/main/results/flagship7e01.png`](cases/main/results/flagship7e01.png) |
| Re = 100   | [`cases/main/results/flagship1e02.png`](cases/main/results/flagship1e02.png) |
| Re = 1000  | [`cases/main/results/flagship1e03.png`](cases/main/results/flagship1e03.png) |
| Re = 10000 | [`cases/main/results/flagship1e04.png`](cases/main/results/flagship1e04.png) |

### Visualising in ParaView

1. **Open** `flagship.pvd` → press Play 
2. **Velocity magnitude**: Filters → Calculator → `mag(velocity)`, colour by result 
3. **Vorticity** (z-component): Filters → Gradient → component (1,0) of ∇u 
4. **Pressure coefficient**: colour by `pressure` field 
5. **Animation**: View → Animation, export as MP4

For physics interpretation, prefer true `velocity` magnitude and/or vorticity fields. `velocity_rgb` is a visualization color field and its magnitude is not a physical speed.

---

## External Resources & References

**Code Repositories & Documentation:**
- [OpenFOAM Tutorials](https://www.openfoam.com/documentation/tutorial-guide)
- [CFD Online Wiki & Forums](https://www.cfd-online.com/)
- [NIST CFD Validation Database](https://www.nist.gov/programs-projects/cfdsupport)
**Visual References:**
- [YouTube: CFD Cavity Flow Visualizations](https://www.youtube.com/results?search_query=CFD+cavity+flow+visualization)
- [YouTube: Cylinder Flow Vortex Shedding](https://www.youtube.com/results?search_query=cylinder+flow+visualization+CFD)
- [OpenFOAM Case Gallery](https://www.openfoam.com/community/case-studies)

**Publication Database:**
- [ResearchGate — CFD Papers](https://www.researchgate.net/)
- [Google Scholar — CFD](https://scholar.google.com/scholar?q=computational+fluid+dynamics)
- [Journal of Computational Physics](https://www.sciencedirect.com/journal/journal-of-computational-physics)

**Key Benchmark Papers:**
- [Ghia et al. (1982) — High-Re solutions for incompressible flow](https://www.researchgate.net/publication/267886860_High-Re_solutions_for_incompressible_flow_using_the_Navier-Stokes_equations_and_a_multigrid_method)
- [Williamson (1996) — Vortex dynamics in the cylinder wake](https://www.annualreviews.org/content/26/1/413)
- [Erturk & Gokcol (2007) — 2-D steady flow in driven skew cavities](https://www.sciencedirect.com/journal/journal-of-computational-physics)
- [Churchill & Chu (1975) — Natural convection correlations](https://www.tandfonline.com/doi/abs/10.1080/15567811.2013.812883)

------------------------------------------------------------------------

## Key Features

**Phase 1 [Complete] Complete (Laminar Flow):**
- [Complete] Staggered grid SIMPLE solver (steady + transient)
- [Complete] Lid-driven cavity benchmark (Re=100, validated vs Ghia)
- [Complete] Energy transport with MMS verification
- [Complete] Async parallel I/O (20-50% speedup)
- [Complete] Animation resume with field loading
- [Complete] Five convection schemes with flux limiters

**Phase 2 [Complete] Complete (Turbulent Flow):**
- [Complete] k-ε turbulence model (production, dissipation, diffusion)
- [Complete] Logarithmic wall functions (high-Re treatment)
- [Complete] Wall y+ diagnostics and validation
- [Complete] High-Re cavity test case (Re=1000+)
- [Complete] Turbulent channel flow (pressure gradient driven)
- [Complete] Source term support in SIMPLE
- [Complete] Comprehensive turbulence documentation

**For Users:**
- Simple CLI interface: `./cfd --animate --name mytest --numTimeSteps 10000`
- ParaView visualization (PVD collections)
- CSV benchmark output (centerline profiles)
- Validation against published data (Ghia et al.)

**For Developers:**
- Clear code architecture (core → numerical → physics → solver)
- Easy to add new cases (template in `cases/`)
- Easy to add new physics (generic scalar transport engine)
- GPU acceleration ready (CUDA backend optional)

---

## Highlights: Animation Resume & Parallel I/O

**Problem:** Creating 10,000 animation frames took minutes due to file I/O bottleneck.

**Solution:**
1. **Async file writing** — Background thread writes while solver computes (20-50% speedup)
2. **Resume capability** — Continue from frame 5000 → 10000 without restart
3. **Field loading** — Resume loads u/v/p state from last VTS file (proper physics)
4. **PVD preservation** — All frames in one PVD file, no metadata duplication

**Usage:**
```bash
# Run 5000 frames
./cfd --animate --name mytest --numTimeSteps 5000

# Later: continue to 10000 (resumes from frame 5000)
./cfd --animate --name mytest --numTimeSteps 10000
# Result: mytest/ has 10000 frames, single mytest.pvd with all entries
```

See: [docs/06-user-guide/animation-guide.md](docs/06-user-guide/animation-guide.md)

------------------------------------------------------------------------

## Governing Equations

**Continuity (Mass conservation):**
$$\nabla \cdot \vec{u} = 0$$

**Momentum (Navier-Stokes):**
$$\frac{\partial \vec{u}}{\partial t} + \vec{u} \cdot \nabla \vec{u} = -\nabla p + \nu \nabla^2 \vec{u}$$

**Generic scalar transport (energy/species):**
$$\frac{\partial (\rho \phi)}{\partial t} + \nabla \cdot (\rho \vec{u} \phi) = \nabla \cdot (\Gamma \nabla \phi) + S_\phi$$

Discretized using finite volume method on **staggered grid** (velocities at faces, pressure at cells).

See: [docs/01-governing-equations/](docs/01-governing-equations/) for full details with derivations.

------------------------------------------------------------------------

## Numerical Method

**Staggered (MAC) grid:**
- Scalar fields at cell centers
- Velocity components at face centers
- Prevents checkerboard pressure oscillations

**Finite Volume discretization:**
- Upwind convection scheme (bounded, O(h) accuracy)
- Central difference diffusion (O(h²) accuracy)
- Five-point Laplacian stencil

**SIMPLE algorithm:**
1. Guess initial pressure
2. Solve momentum equations → intermediate velocity
3. Solve pressure-correction equation
4. Correct pressure and velocity
5. Iterate until convergence (typically 5-10 iterations per timestep)

See: [docs/02-numerical-methods/](docs/02-numerical-methods/) for discretization details.

------------------------------------------------------------------------

## Build Instructions

**Requirements:**
- CMake ≥ 3.16
- C++17 compiler (g++, clang, MSVC)
- Eigen ≥ 3.3 (linear algebra)
- CUDA ≥ 11 (optional, for GPU acceleration)
- BLAS/LAPACK (optional, for sparse solvers)

**Build:**
```bash
mkdir build && cd build
cmake .. && make -j4
```

**With GPU (CUDA):**
```bash
cmake -DCFD_BACKEND=CUDA .. && make -j4
```

**Run tests:**
```bash
./cfd --help # Show options
./cfd # Steady-state lid cavity (Re=100)
./cfd --animate --numTimeSteps 50 # Transient animation (50 frames)
./cfd --case energy_transport # Validation test
```

**Output files:**
- `animation/` folder (if `--animate`)
- `animation.pvd` (ParaView collection file)
- `lid_driven_cavity.vtk` (final solution)
- `centerline_u.csv`, `centerline_v.csv` (benchmarks)

See: [docs/06-user-guide/running-simulations.md](docs/06-user-guide/running-simulations.md) for detailed setup.

------------------------------------------------------------------------

## Validation & Benchmarks

**Lid-Driven Cavity (Re=100):**
- Compared against Ghia et al. (1982) benchmark
- Centerline L2 error: 0.344 (u), 0.226 (v)
- Mass residual: ~10⁻¹⁹
- Status: [Complete] **PASS**

**Energy Transport (Manufactured Solution):**
- Analytical solution known exactly
- L2 error: 0.00597 (0.60% relative)
- Spatial convergence: O(h²) verified
- Status: [Complete] **PASS**

See: [docs/05-cases-validation/](docs/05-cases-validation/) for full validation methodology and results.

------------------------------------------------------------------------

## Project Structure

```
cfd-engine/
├── core/ → Mesh, Field, UField, VField, BoundaryCondition
├── numerics/ → ConvectionSchemes, Laplacian, Gradient, Divergence
├── assembly/ → Assemblers for momentum, pressure, scalars
├── physics/ → NavierStokes, ScalarTransport, TurbulenceKEpsilon
├── solver/ → SIMPLE, StaggeredSIMPLE, TimeLoop
├── io/ → VTKWriter, PVDWriter (async+resume!), VTSReader
├── cases/ → LidDrivenCavity, EnergyTransport, test cases
├── config/ → SolverConfig (all parameters)
├── linalg/ → Linear solver backends (Eigen, CUDA)
├── validation/ → Ghia reference data
├── docs/ → **NEW: Comprehensive documentation (7 sections)**
│ ├── 00-course-foundation/
│ ├── 01-governing-equations/
│ ├── 02-numerical-methods/
│ ├── 03-physics-modules/
│ ├── 04-implementation/
│ ├── 05-cases-validation/
│ ├── 06-user-guide/
│ ├── 07-future-directions/
│ └── INDEX.md
├── slides/ → Original course materials (11 PDF files)
├── README.md → **You are here**
├── ARCHITECTURE.md → Design and layer documentation
├── CHANGELOG.md → Version history and changes
└── CMakeLists.txt

See: [docs/04-implementation/project-structure.md](docs/04-implementation/project-structure.md) for detailed explanation of each module.
```

------------------------------------------------------------------------

## Verification Strategy

Phase 1 validation covers:

- **Lid-driven cavity** against Ghia et al. (1982) benchmark
- **Energy transport** manufactured solution (Method of Manufactured Solutions — MMS)
- **Residual and divergence** norm monitoring
- Grid refinement studies (future phases)

See: [docs/05-cases-validation/convergence-studies.md](docs/05-cases-validation/convergence-studies.md) for detailed explanation.

------------------------------------------------------------------------

## Design Principles

This solver is designed to:

- **Make every numerical operator explicit** — See exactly what equation is discretized
- **Avoid hidden abstractions** — Easy to debug and understand
- **Support modular extension** — Add new physics/cases without refactoring
- **Provide full transparency** — Clear indexing, no magic
- **Maintain educational value** — Learn CFD through working code

The code prioritizes **clarity and extensibility** over industrial complexity.

------------------------------------------------------------------------

## Roadmap

**Phase 1 [Complete] (Complete)**
- SIMPLE solver + validation + async I/O + resume
- Comprehensive documentation aligned with course lectures

**Phase 2 Plan (Planned)**
- Turbulence modeling (k-ε, LES)
- Grid independence studies
- More validation cases

**Phase 3+ (Vision)**
- Multiphase flows
- 3D support
- Advanced preconditioners
- Python API

See: [docs/07-future-directions/roadmap.md](docs/07-future-directions/roadmap.md) for details.

------------------------------------------------------------------------

## Learning Resources

**Lectures (from course):**
- Intro_marked.pdf — CFD fundamentals
- GovEqI/II_marked.pdf — Governing equations
- FVII_unmarked.pdf — Finite Volume discretization
- SpecHeat_marked.pdf — Species & heat transport
- TurbI/II/III.pdf — Turbulence modeling
- MultiPhase1/2.pdf — Multiphase flows (future)

Located in: `/slides/` folder

**Code + Docs:**
- Run: `./cfd --help` (see all CLI options)
- Explore: `cases/` folder (reference implementations)
- Read: `docs/INDEX.md` (comprehensive 50+-page documentation)
- Start with: [docs/00-course-foundation/learning-path.md](docs/00-course-foundation/learning-path.md)

------------------------------------------------------------------------

## Contributing

Want to extend the solver? See:
- [docs/04-implementation/adding-cases.md](docs/04-implementation/adding-cases.md) — Add new test case
- [docs/04-implementation/adding-physics.md](docs/04-implementation/adding-physics.md) — Add new physics module
- [docs/07-future-directions/](docs/07-future-directions/) — Feature backlog

------------------------------------------------------------------------

## License

TBD (Educational use assumed)

------------------------------------------------------------------------

## Next Steps

→ **New to CFD?** Start with [docs/00-course-foundation/overview.md](docs/00-course-foundation/overview.md)

→ **Ready to run?** Go to [docs/06-user-guide/running-simulations.md](docs/06-user-guide/running-simulations.md)

→ **Want to learn deeply?** Follow [docs/00-course-foundation/learning-path.md](docs/00-course-foundation/learning-path.md) (50 hours, 9 levels)

→ **See everything:** [docs/INDEX.md](docs/INDEX.md)