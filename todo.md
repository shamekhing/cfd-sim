# CFD Engine - Phase 3 Development Roadmap

**Last Updated:** March 7, 2026 
**Current Status:** Phase 2 Complete (k-ε turbulence + CUDA backend) 
**Next Phase:** Advanced Solvers & Physics Extensions

---

## Plan Executive Summary

This document outlines the progression from Phase 2 (current) to Phase 3, with three strategic paths:

| Goal | Path | Duration | Complexity |
|------|------|----------|-----------|
| Learning Course completion | Validation + Multigrid | 4 weeks | Medium |
| Research Research/Thesis | LES or k-ω SST | 5-8 weeks | High |
| Performance Production ready | Multiphase (VOF) + GPU opt. | 8-12 weeks | Very High |
| Lecture Learning focused | Validation → Multigrid → LES | 8-10 weeks | Medium-High |

---

## Focus Your Options at a Glance

```
 VALIDATION PASS (Phase 3.0)
 ↓
 ┌─────────────────┼─────────────────┐
 ↓ ↓ ↓
 MULTIGRID GPU OPT. [SKIP]
 Solver & Profiling
 ↓ ↓
 └─────────────────┼─────────────────┘
 ↓
 PICK ONE PHASE 3 EXTENSION
 ┌─────────────────┼─────────────────┐
 ↓ ↓ ↓
 LES k-ω SST VOF
 (Best for (Industry (Novel
 Learning) Standard) Physics)
```

---

## Acceleration PHASE 3.0: Validation & Benchmarking (Weeks 1-2)

**Status:** Essential first step | **Effort:** 40-50 hours | **Priority:** High PriorityHigh PriorityHigh Priority

### Objectives
- [Complete] Verify k-ε solver against published benchmarks
- [Complete] Perform systematic mesh convergence study
- [Complete] Quantify discretization error
- [Complete] Document validation report

### Detailed Tasks

#### 1.1 Mesh Convergence Study
**Files involved:** `cases/examples/LidDrivenCavity.hpp`, `docs/TECHNICAL_REFERENCE.md`

- Run 5 mesh refinements: 32×32 → 64×64 → 128×128 → 256×256 → 512×512
- Track variables: $u_{max}$, $v_{max}$, center pressures
- Compute error vs finest mesh
- Plot convergence rate (should be $O(\Delta x^2)$ for 2nd order)
- **Deliverable:** Convergence table + plot in `docs/validation/convergence.md`

**Code changes needed:**
```cpp
// Add to SolverConfig
struct ValidationConfig {
 std::vector<int> mesh_levels = {32, 64, 128, 256};
 int num_runs = 4;
 bool output_convergence_report = true;
};
```

#### 1.2 Lid-Driven Cavity Validation (Re=100, 1000)
**Reference:** Ghia et al. (1982) benchmark

- Compare $u(0.5, y)$ and $v(x, 0.5)$ with Ghia data
- Generate L2 error norm
- Document deviations
- **Deliverable:** `docs/validation/ghia_comparison.csv` + plot

**Implementation location:** `src/validation/GhiaComparison.cpp` (new)

#### 1.3 Cylinder Vortex Shedding Validation
**Reference:** Williamson (1996)

- Validate Strouhal number at Re=100: $St_{Williamson} = 0.164 \approx St_{computed}$
- Validate drag coefficient: $C_d \approx 1.3$
- Run 3 Re values: 100, 150, 200
- **Deliverable:** `docs/validation/cylinder_validation.md`

#### 1.4 Energy Transport MMS Convergence
**Already implemented** - just extend documentation

- Verify $L_2$ error scales as $O(\Delta x^2)$
- Add thermal error analysis to `docs/TECHNICAL_REFERENCE.md`

#### 1.5 Turbulent Channel Flow
**Reference:** DNS database (Moser et al. 1999)

- Compare mean velocity profile $\bar{u}^+$ vs DNS
- Compare Reynolds stress $\overline{u'v'}$ vs DNS
- **Deliverable:** `docs/validation/channel_dns_comparison.md`

---

## PHASE 3.1: Advanced Linear Solvers (Weeks 2-4)

**Status:** High-impact performance upgrade | **Effort:** 50-60 hours | **Priority:** High PriorityHigh PriorityHigh Priority

### Objectives
- [Complete] Implement Multigrid V-cycle (10-50× speedup on Poisson)
- [Complete] Add GMRES for general systems
- [Complete] Benchmark solver performance
- [Complete] Integrate with SIMPLE

### Detailed Tasks

#### 3.1.1 Multigrid Implementation (V-cycle)
**New files:**
- `include/linear/Multigrid.hpp`
- `src/linear/Multigrid.cpp`
- `include/linear/Restriction.hpp`
- `include/linear/Prolongation.hpp`

**Implementation outline:**
```
Multigrid V-cycle for Ax = b:
1. Pre-smooth: x ← x + ω(b - Ax) [2-3 iterations]
2. Restrict: r_coarse = P^T(b - Ax)
3. Solve coarse: e_coarse = AMG(r_coarse)
4. Prolongate: x ← x + P·e_coarse
5. Post-smooth: x ← x + ω(b - Ax) [2-3 iterations]
6. Repeat 1-5 for <5 V-cycles or until converged
```

**Key classes to modify:**
- `include/linear/ILinearSolver.hpp` - add `setCoarseSolver()` interface
- `solver/StaggeredSIMPLE.hpp` - replace CG/BiCGSTAB with Multigrid option

**Success criteria:**
- 10-50× speedup on pressure Poisson
- Residual reduction similar to CG
- Scales well to 1000×1000 grids

#### 3.1.2 GMRES Implementation (Restarted GMRES-m)
**New file:** `include/linear/GMRES.hpp`

- For general non-symmetric matrices
- Use 30 restart length
- Integrate with momentum solvers

#### 3.1.3 ILU Preconditioner
**New file:** `include/linear/ILUPreconditioner.hpp`

- Incomplete LU factorization
- Use as preconditioner for CG/GMRES
- Reduce iteration count by 30-50%

#### 3.1.4 Benchmarking & Integration
**New file:** `tests/validation/SolverBenchmark.cpp`

- Compare CG vs Multigrid on 64×64 to 512×512 grids
- Measure wall time, iterations, memory
- Create solver selection guide
- **Deliverable:** `docs/solver-performance-report.md`

---

## Metrics PHASE 3.2: GPU Memory Optimization (Weeks 3-5)

**Status:** Optional but recommended | **Effort:** 30-40 hours | **Priority:** High PriorityHigh Priority

### Objectives
- [Complete] Profile CUDA kernels (identify bottlenecks)
- [Complete] Implement memory pooling
- [Complete] Optimize kernel launch parameters
- [Complete] Achieve 2-5× speedup on large grids

### Detailed Tasks

#### 3.2.1 CUDA Profiling
**Tools:** `nsys`, `ncu` (NVIDIA profilers)

- Profile assembly kernels with 256×256 grid
- Identify memory bandwidth bottlenecks
- Check occupancy, warp efficiency
- **Deliverable:** `docs/gpu-profile-baseline.md`

#### 3.2.2 Memory Pooling
**File:** `backend/cuda/include/CudaMemoryPool.hpp` (new)

- Pre-allocate memory pools at startup
- Reuse allocations to avoid malloc/free overhead
- Reduce host-device sync points

#### 3.2.3 Kernel Optimization
**Files:** `backend/cuda/src/CudaAssembly.cu`, `backend/cuda/src/CudaLaplacian.cu`

- Tune block/thread sizing
- Minimize register pressure
- Use shared memory for stencils
- **Expected gain:** 2-5× on assembly operations

---

## Flow PHASE 3.3: Advanced Physics (Pick ONE)

Choose based on your goal and interests:

---

### **OPTION A: LES (Large Eddy Simulation)** Learning

**Status:** Best for learning | **Effort:** 60-80 hours | **Priority:** High PriorityHigh PriorityHigh Priority

#### Why LES?
- [Complete] Simpler than k-ω SST (fewer equations)
- [Complete] Captures large eddy dynamics
- [Complete] Great physics learning
- [Complete] Production use in industry

#### Implementation Roadmap

**Phase 3.3.A1: Smagorinsky Model (Weeks 4-6)**

**New files:**
- `include/physics/LES.hpp`
- `src/physics/LES.cpp`
- `include/numerics/FilterOperator.hpp`

**Theory:**
$$\nu_t = (C_s \Delta)^2 |\bar{S}|$$

where:
- $\Delta = (\Delta x \Delta y)^{1/2}$ (filter width)
- $|\bar{S}| = \sqrt{2\bar{S}_{ij}\bar{S}_{ij}}$ (strain rate)
- $C_s \approx 0.1-0.2$ (Smagorinsky constant)

**Implementation steps:**
1. Implement test filter (Gaussian or box)
2. Compute subgrid scale viscosity
3. Assemble momentum equation with SGS stress
4. Validate on homogeneous isotropic turbulence (HIT)

**Code location:**
```cpp
namespace physics {
 class LES {
 public:
 void initialize(const core::Mesh& mesh);
 void computeSubgridViscosity(
 const core::UField& u,
 const core::VField& v,
 const core::Mesh& mesh);
 core::Field& nu_t() { return nu_t_; }
 private:
 core::Field nu_t_;
 linear::Scalar C_s_ = 0.1;
 };
}
```

**Test cases:**
- Homogeneous isotropic turbulence (HIT)
- Channel flow at Re_τ=5600 (compare to DNS)
- High-Re lid-driven cavity (Re=10,000)

**Deliverables:**
- `include/physics/LES.hpp` - Smagorinsky implementation
- `cases/examples/HighReynoldsCavityLES.hpp` - LES case
- `docs/les-implementation.md` - theoretical background
- Validation report vs DNS/experiments

**Success criteria:**
- Energy spectrum matches DNS in inertial range
- Mean velocity profiles within 10% of DNS
- Computational cost < k-ε + improved turbulence capture

---

### **OPTION B: k-ω SST Turbulence Model** Focus

**Status:** Industry standard | **Effort:** 50-70 hours | **Priority:** High PriorityHigh PriorityHigh Priority

#### Why k-ω SST?
- [Complete] Better near-wall than k-ε (no wall functions needed)
- [Complete] Industry standard (ANSYS, OpenFOAM)
- [Complete] Direct k-ε replacement
- [Complete] Blends k-ε (far-field) + k-ω (near-wall)

#### Implementation Roadmap

**Phase 3.3.B1: k-ω SST Model (Weeks 4-6)**

**New files:**
- `include/physics/TurbulenceKOmegaSST.hpp`
- `src/physics/TurbulenceKOmegaSST.cpp`
- `include/config/KOmegaSSTConfig.hpp`

**Theory:**
$$\frac{\partial k}{\partial t} + u_j \nabla_j k = P_k - \beta^* \rho \omega k + \nabla_j[(\nu + \sigma_k \nu_t)\nabla_j k]$$

$$\frac{\partial \omega}{\partial t} + u_j \nabla_j \omega = \frac{\gamma}{\nu_t}P_k - \beta \rho \omega^2 + \nabla_j[(\nu + \sigma_\omega \nu_t)\nabla_j \omega] + 2(1-F_1)\frac{\sigma_{\omega,2}}{\omega}\nabla_k \cdot \nabla\omega$$

Where:
- $F_1, F_2$ = blending functions (switch k-ε ↔ k-ω near walls)
- $\nu_t = a_1 k / \max(a_1 \omega, \Omega F_2)$
- Constants: $\beta^*=0.09$, $\beta_0=0.075$, etc.

**Implementation steps:**
1. Blend k-ε (far-field) + k-ω (near-wall)
2. Implement blending functions $F_1, F_2$
3. Solve without wall functions (resolve to $y^+ \sim 0.1$)
4. Validation on boundary layers

**Test cases:**
- Flat plate turbulent boundary layer
- Channel flow with low-Re resolution (y+=0.1)
- Cylinder wake Re=140 (compare to experiments)
- High-Re cavity (Re=10,000)

**Deliverables:**
- `include/physics/TurbulenceKOmegaSST.hpp` - full implementation
- `include/config/KOmegaSSTConfig.hpp` - model constants
- `cases/examples/FlatPlateBoundaryLayer.hpp` - validation case
- `docs/k-omega-sst-implementation.md` - theory + implementation
- Performance report vs k-ε

**Success criteria:**
- No wall functions needed
- Skin friction coefficient within 5% of experiments
- Better near-wall accuracy than k-ε
- Handles adverse pressure gradients

---

### **OPTION C: VOF Multiphase (Free Surface)** Multiphase

**Status:** Most novel | **Effort:** 80-100 hours | **Priority:** High PriorityHigh Priority

#### Why VOF?
- [Complete] Most complex (great for thesis)
- [Complete] Real-world applications (dam break, sloshing)
- [Complete] Combines new physics (advection, surface tension)
- [Complete] High publication potential

#### Implementation Roadmap

**Phase 3.3.C1: VOF Advection (Weeks 4-6)**

**New files:**
- `include/physics/VOF.hpp`
- `src/physics/VOF.cpp`
- `include/numerics/InterfaceTracking.hpp`
- `include/numerics/SurfaceTension.hpp`

**Theory:**
$$\frac{\partial f}{\partial t} + u_j \nabla_j f = 0$$

where $f$ = volume fraction (1 = liquid, 0 = gas)

**Implementation steps:**
1. PLIC (Piecewise Linear Interface Calculation)
2. Flux-corrected transport (FCT) for advection
3. Surface tension modeling (CSF - Continuum Surface Force)
4. Adaptive mesh refinement near interface

**Test cases:**
- Rayleigh-Taylor instability
- Rising bubble
- Dam break (experimental validation)
- Sloshing tank

**Deliverables:**
- `include/physics/VOF.hpp` - advection solver
- `include/numerics/SurfaceTension.hpp` - surface tension
- `cases/examples/RisingBubble.hpp` - validation case
- `docs/vof-implementation.md` - complete writeup
- Experimental comparison video

**Success criteria:**
- Mass conservation error < 1%
- Interface stays sharp (no diffusion)
- Rising bubble terminal velocity matches theory
- Good performance on fine grids

---

## Timeline Detailed Timeline Recommendations

| Path | Week 1-2 | Week 3-4 | Week 5-6 | Week 7-8 |
|------|----------|----------|----------|----------|
| **Learning** | Validation | Multigrid | LES (Phase 1) | LES (Phase 2) |
| **Research** | Validation | GPU Opt. | k-ω SST | k-ω SST refinement |
| **Industrial** | Validation | Multigrid | GPU Opt. | VOF (Phase 1) |
| **Quick Win** | Validation | Multigrid | — | Deploy & benchmark |

---

## [Complete] Success Criteria Matrix

| Milestone | Validation | Multigrid | LES | k-ω SST | VOF |
|-----------|-----------|-----------|-----|---------|-----|
| Code complete | | Y/N | | | |
| Tests pass | | 2-5× faster | Energy spectrum | Near-wall | Mass cons. |
| Documented | | Benchmark report | Theory + impl. | Theory + impl. | Theory + videos |
| Published | Optional | Optional | High | Medium | High |

---

## Implementation Implementation Checklist Template

```markdown
# [Feature Name] Implementation Checklist

## Phase 1: Foundation
- [ ] Create header files + stubs
- [ ] Write unit tests (stub implementations)
- [ ] Set up CI/CD
- [ ] Document interfaces

## Phase 2: Core Implementation
- [ ] Implement core algorithms
- [ ] Pass unit tests
- [ ] Integrate with SIMPLE solver
- [ ] Add convergence monitoring

## Phase 3: Validation
- [ ] Run benchmark cases
- [ ] Compare vs reference data
- [ ] Document results
- [ ] Performance profiling

## Phase 4: Production
- [ ] Code review
- [ ] Add regression tests
- [ ] Optimize hotspots
- [ ] Write user guide
```

---

## Learning Key References by Topic

### Validation & Benchmarking
- **Ghia et al. (1982)** — Lid-driven cavity reference solution
- **Williamson (1996)** — Vortex shedding, Strouhal number
- **Moser et al. (1999)** — Channel flow DNS database

### Advanced Solvers
- **Multigrid:** Briggs, Henson, McCormick (2000) "A Multigrid Tutorial"
- **GMRES:** Saad & Schultz (1986)
- **ILU Preconditioner:** Benzi (2002)

### Advanced Physics
- **LES:** Sagaut (2006) "Large Eddy Simulation"
- **k-ω SST:** Menter (1994) "Two-equation eddy-viscosity turbulence models"
- **VOF:** Hirt & Nichols (1981), Scardovelli & Zaleski (1999)

---

## Focus Quick Decision Tree

```
START
 ↓
Have you validated Phase 2? → NO → Do PHASE 3.0 first
 ↓ YES
 
What's your PRIMARY goal?
 ├─→ Learn CFD fundamentals? → VALIDATION + MULTIGRID
 ├─→ Build thesis/publication? → LES or k-ω SST
 ├─→ Get production-ready? → MULTIGRID + k-ω SST
 └─→ Explore new physics? → VOF + Multiphase

What's your time budget?
 ├─→ 4 weeks? → Validation + Multigrid
 ├─→ 6-8 weeks? → + One advanced physics
 └─→ 12+ weeks? → Full Phase 3 + GPU opt + advanced physics
```

---

## Support Support & Resources

- **Code Review:** Check existing implementations in `include/physics/`
- **Documentation:** `docs/TECHNICAL_REFERENCE.md`
- **Testing:** Run with `./run_all_cases.sh`
- **Benchmarks:** Output in `build_cuda/output/`

---

## Last Updated
March 7, 2026 | Next review: After Phase 3.0 completion