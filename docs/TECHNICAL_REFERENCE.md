# Technical Reference: Code Architecture & Implementation

Complete implementation details for developers extending the CFD solver.

**Quick Navigation:**
- [Architecture Overview](#architecture-overview) — 4-layer design
- [Core Infrastructure](#core-infrastructure) — Mesh, fields, indexing 
- [Staggered Grid](#staggered-grid) — MAC scheme
- [Discretization](#discretization-layer) — Operators and schemes
- [SIMPLE Algorithm](#simple-algorithm) — Pressure-velocity coupling
- [Turbulence](#turbulence-implementation) — k-ε model and wall functions
- [Assembly](#assembly-layer) — Matrix construction
- [Innovations](#solver-innovations) — Quasi-transient, direct SparseLU

**See also:**
- Learning & theory → [COURSE.md](COURSE.md)
- Building & running → [USER_GUIDE.md](USER_GUIDE.md) 
- Project overview → [README.md](README.md)

---

## Architecture Overview

**4-layer design:**

```
┌────────────────────────────────────────────┐
│ Physics (SIMPLE, Turbulence) │
├────────────────────────────────────────────┤
│ Assembly (Matrix construction) │
├────────────────────────────────────────────┤
│ Discretization (Operators, schemes) │
├────────────────────────────────────────────┤
│ Core (Mesh, fields, indexing) │
└────────────────────────────────────────────┘
```

### Layer 1: Core Infrastructure

**Components:**
- `core/Mesh.hpp` — Cartesian structured grid
- `core/Field.hpp` — Cell-centered scalars
- `core/UField.hpp`, `core/VField.hpp` — Staggered velocities
- `linalg/` — Sparse matrices, solvers

### Layer 2: Discretization 

**Components:**
- `numerics/Laplacian.hpp` — 5-point Laplacian
- `numerics/Gradient.hpp` — Pressure gradient
- `numerics/Divergence.hpp` — Velocity divergence
- `numerics/ConvectionSchemes.cpp` — 5 schemes + limiters

### Layer 3: Assembly

**Components:**
- `assembly/StaggeredMomentumAssembler.hpp` — Momentum equations
- `assembly/StaggeredPressureCorrectionAssembler.hpp` — Pressure Poisson
- `assembly/GenericScalarAssembler.hpp` — Convection-diffusion

### Layer 4: Physics

**Components:**
- `solver/StaggeredSIMPLE.hpp` — Pressure-velocity coupling
- `physics/TurbulenceKEpsilon.hpp` — k-ε turbulence model
- `physics/ScalarTransport.hpp` — Generic scalars

---

## Core Infrastructure

### Mesh

**Uniform Cartesian grid:**
```cpp
core::Mesh mesh(Nx, Ny, Lx, Ly);
double dx = Lx / Nx;
double dy = Ly / Ny;
```

**Cell indices:** $\text{id} = i + j \cdot N_x$

### Staggered Grid

**Locations:**
- **Pressure (p):** Cell centers $(i, j)$
- **U-velocity (u):** Vertical faces $(i+0.5, j)$
- **V-velocity (v):** Horizontal faces $(i, j+0.5)$

**Why:** Prevents checkerboard pressure oscillations, natural divergence

```
+-------+-------+
│ p│u │ p │
│ v u │ v u │ p = cell center
│────────────────│ u = vertical face
│ p│u │ p │ v = horizontal face
+-------+-------+
```

### Indexing

```cpp
// Cell access
Field<double> p(mesh);
p(i, j) = value; // Stored as p[i + j*Nx]

// Face velocities
UField u(mesh); // u at (i+0.5, j)
VField v(mesh); // v at (i, j+0.5)
```

---

## Discretization Layer

### Laplacian (5-point stencil)

$$\nabla^2 \phi_P = \frac{\phi_E - 2\phi_P + \phi_W}{\Delta x^2} + \frac{\phi_N - 2\phi_P + \phi_S}{\Delta y^2}$$

**Implementation:** `numerics/Laplacian.hpp`

### Convection Schemes

**5 options in `numerics/ConvectionSchemes.cpp`:**

| Scheme | Order | Stability | Use |
|--------|-------|-----------|-----|
| Upwind | 1 | Robust | Default (bounded) |
| Central | 2 | Conditional | Low-Re only |
| QUICK | 3 | Limited | Smooth regions |
| MUSCL | 2 TVD | TVD | Shock capture |
| Van Leer | 2 | Smooth | General |

**Default:** Bounded upwind (ensures diagonal dominance)

### Bounded Upwind

**Ensures:** $a_P \geq |a_E| + |a_W| + |a_N| + |a_S|$

```cpp
double uC = 0.5 * (u(i, j) + u(i+1, j));
if (uC >= 0.0) {
 aP += uC / dx;
 aW += uC / dx;
} else {
 aP += (-uC) / dx;
 aE += (-uC) / dx;
}
```

---

## SIMPLE Algorithm

**Pressure-Velocity Coupling:**

```
1. Guess p*
2. Solve momentum → u*, v* (violates mass conservation)
3. Compute divergence: div(u*)
4. Solve pressure correction: ∇·(1/aP·∇p') = div(u*)
5. Correct: p_new = p* + αp·p'
 u_new = u* - correction(∇p')
6. Iterate until convergence
```

**Key aspects:**
- Quasi-transient stabilization: Add $\Delta t_{\text{pseudo}} = 0.1 \Delta x$ term
- Direct SparseLU solver (not BiCGSTAB)
- Bounded upwind for diagonal dominance
- Underrelaxation: $\alpha_p \approx 0.2-0.3$, $\alpha_u \approx 0.5-0.7$

**Implementation:** `solver/StaggeredSIMPLE.hpp`

---

## Turbulence Implementation

### k-ε Model

**Turbulent kinetic energy:**
$$\frac{\partial k}{\partial t} + \vec{u} \cdot \nabla k = \nabla \cdot[(\nu + \nu_t/\sigma_k) \nabla k] + P_k - \varepsilon$$

**Dissipation:**
$$\frac{\partial \varepsilon}{\partial t} + \vec{u} \cdot \nabla \varepsilon = \nabla \cdot[(\nu + \nu_t/\sigma_\varepsilon) \nabla \varepsilon] + C_{1\varepsilon}\frac{\varepsilon}{k}P_k - C_{2\varepsilon}\frac{\varepsilon^2}{k}$$

**Turbulent viscosity:**
$$\nu_t = C_\mu \frac{k^2}{\varepsilon}$$

**Production:**
$$P_k = \nu_t \cdot 2\|\mathbf{S}\|^2$$

**Constants:** $C_\mu=0.09$, $\sigma_k=1.0$, $\sigma_\varepsilon=1.3$, $C_{1\varepsilon}=1.44$, $C_{2\varepsilon}=1.92$

### Wall Functions

**Logarithmic law (y+ > 11.25):**
$$u^+ = \frac{1}{\kappa}\ln(E \cdot y^+)$$

Where: $u^+ = u/u_\tau$, $y^+ = y \cdot u_\tau/\nu$, $\kappa=0.41$, $E=9.8$

**Newton iteration** to solve for $u_\tau$ implicitly (3-5 iterations)

**Wall boundary conditions:**
- If $y^+ > 11.25$: $k_w = u_\tau^2/\sqrt{C_\mu}$, $\varepsilon_w = u_\tau^3/(\kappa y_p)$
- Else: $k_w = 0$, $\varepsilon_w = 2\nu k_p / y_p^2$

**File:** `physics/TurbulenceKEpsilon.hpp::applyWallFunctions()`

### y+ Diagnostics

Compute dimensionless wall distance for all walls:
$$y^+ = \frac{y \cdot u_\tau}{\nu}$$

**File:** `solver/ResidualMonitor.hpp::wallYPlus()`

**Guidelines:**
- y+ < 5: Viscous sublayer (wall-resolved)
- 5 < y+ < 30: Buffer layer (wall functions less accurate)
- y+ > 30: Log layer (wall functions fully valid)

---

## Assembly Layer

### Momentum Assembly

**U-momentum equation at face (i, j):**
$$a_P u_P = \sum a_{\text{nb}} u_{\text{nb}} + b - \frac{\Delta y}{\Delta x}(p_E - p_P) + S_u$$

**File:** `assembly/StaggeredMomentumAssembler.hpp`

**Returns:** Linear system $(A_u, b_u)$ for $A_u u = b_u$

### Pressure Correction Assembly

**Poisson for pressure correction:**
$$\frac{1}{a_P^u} \frac{p'_E - 2p'_P + p'_W}{\Delta x^2} + \frac{1}{a_P^v} \frac{p'_N - 2p'_P + p'_S}{\Delta y^2} = \text{div}(u^*)$$

**File:** `assembly/StaggeredPressureCorrectionAssembler.hpp`

**Boundary conditions:** Neumann at walls ($\partial p'/\partial n = 0$)

### Generic Scalar Assembly

**Convection-diffusion:**
$$\frac{\partial \phi}{\partial t} + \nabla \cdot (\vec{u}\phi) = \nabla \cdot (\Gamma \nabla \phi) + S_\phi$$

**Features:**
- Implicit Euler (1st-order time)
- Selectable convection scheme
- Central diffusion (2nd-order)

**File:** `assembly/GenericScalarAssembler.hpp`

**Used for:** Temperature (energy), k (turbulence), ε (dissipation), custom scalars

---

## Solver Innovations

### 1. Quasi-Transient SIMPLE

**Problem:** Pure steady-state creates singular momentum matrices at ~160 iterations

**Solution:** Add fictitious timestep $\Delta t_{\text{pseudo}} = 0.1 \Delta x$

**Effect:** Adds diagonal stabilization term $\rho u / \Delta t_{\text{pseudo}}$

**Result:** Robust convergence in 300 iterations without affecting steady-state

### 2. Direct SparseLU

**Problem:** BiCGSTAB iterative solver fails on ill-conditioned staggered systems

**Solution:** Use Eigen::SparseLU (dense factorization)

**Trade:** Deterministic + always succeeds, but O(n³) (not scalable to large grids)

**Suitable:** Research/validation grids (64×64, 128×128)

### 3. Bounded Upwind Momentum

**Ensures diagonal dominance:**
```cpp
// Bounded upwind for momentum
// Guarantees: a_P >= sum(|off-diagonal|)
```

---

## Implementation Status

**Phase 1 (v1.0.0):**
- [Complete] Staggered SIMPLE with quasi-transient
- [Complete] Generic scalar transport
- [Complete] 5 convection schemes
- [Complete] Lid-driven cavity validation
- [Complete] Energy transport MMS test

**Phase 2 (v2.0.0):**
- [Complete] k-ε turbulence model
- [Complete] Wall functions + y+ diagnostics
- [Complete] Source terms in momentum
- [Complete] High-Re cavity validation
- [Complete] Channel flow validation

---

## Extension Points

### Adding Physics

1. Create config: `config/MyPhysicConfig.hpp`
2. Implement: `physics/MyPhysics.hpp`
3. Solve in SIMPLE loop (pass physics pointer)

### New Convection Scheme

1. Extend `ConvectionScheme` enum
2. Add interpolation function
3. Update switch statement in `interpolate()`

### New Turbulence Model (e.g., k-ω SST)

1. Create `physics/TurbulenceKOmegaSST.hpp`
2. Implement: `computeNut()`, `solve()`, `applyWallFunctions()`
3. Use same interface as k-ε

---

## Input/Output System

### Overview

The IO system provides optimized output with error resilience:

**Components:**
- `io/PVDWriter.hpp` — Async VTS frame writing + PVD time-series (for ParaView animation)
- `io/VTSReader.hpp` — Restart support (reads previous VTS file)
- `io/VTKWriter.hpp` — Legacy single-field VTK format
- `io/CSVWriter.hpp` — Convergence and profile export
- `io/BenchmarkOutput.hpp` — Centerline extractions for validation
- `io/Logger.hpp` — Logging utilities (currently minimal; mostly stdout)

### Performance Characteristics

**PVDWriter Optimizations:**
1. **Async Background Writing** — Field data written in separate thread while solver continues; typically 5-10 ms overhead per frame on modern SSDs
2. **Move Semantics** — Field copies moved into async lambda to avoid extra allocations (1× copy instead of 2×)
3. **Efficient Parsing** — VTSReader uses lightweight streaming parser (no regex) for restart operations; ~10× faster than regex-based approach for large grids (>10k points)
4. **Frame Numbering** — 7-digit padding (e.g., `0000000`, `9999999`) supports up to 10 million simulation steps

**I/O Error Policy:**
- File write failures (disk full, permission denied, invalid path) emit `[IO WARNING]` to stderr
- **Simulation continues** despite IO failures (robust error handling)
- CSV writes throw exceptions for validation test compatibility
- Restart failures (unreadable VTS) log warnings; solver starts from initial conditions

### File Format Details

**VTS (XML StructuredGrid):**
- Point coordinates (Float32): x, y, 0
- Pressure scalar field (Float32)
- Velocity magnitude scalar (Float32)
- Velocity vector field (3×Float32)
- Extensible for additional fields (temperature, turbulence K/ε, etc.)
- ASCII format preserves readability and portability

**PVD (Time-series metadata):**
- XML collection linking VTS frames with timestep/time values
- Relative paths for portability
- Loadable directly in ParaView for animation playback

**CSV (Convergence & profiles):**
- Comma-separated values, scientific notation (precision 8 by default)
- Residuals: iteration, mass_residual, u_residual, v_residual, k_residual, eps_residual, y_plus
- Profiles: coordinate (y or x) and field value (u or v)

### Restart Mechanism

1. **Detection** — `PVDWriter::detectExistingFrames()` scans output directory for highest-numbered VTS file
2. **Loading** — If resume requested: call `VTSReader::read(filename, mesh, u, v, p)` to restore field state
3. **Continuation** — Next write starts from frame N+1 (no overwrites)
4. **Robustness** — If restart file unreadable, solver logs warning and starts from initial conditions (continues rather than crashes)

---

## File References

**Core infrastructure:**
- `core/Mesh.hpp`, `core/Field.hpp`, `core/UField.hpp`, `core/VField.hpp`

**Discretization:**
- `numerics/Laplacian.hpp`, `numerics/ConvectionSchemes.cpp`, `numerics/Gradient.hpp`

**Assembly:**
- `assembly/StaggeredMomentumAssembler.hpp`, `assembly/StaggeredPressureCorrectionAssembler.hpp`, `assembly/GenericScalarAssembler.hpp`

**Solvers:**
- `solver/StaggeredSIMPLE.hpp`, `physics/TurbulenceKEpsilon.hpp`, `physics/ScalarTransport.hpp`

---

## See Also

- **[COURSE.md](COURSE.md)** — Theory and learning path
- **[USER_GUIDE.md](USER_GUIDE.md)** — Building, running, validation
- **[README.md](README.md)** — Quick start
