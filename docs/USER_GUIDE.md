# User Guide: Building, Running & Validation

Complete guide for building the solver, running simulations, and validating results.

**Quick Start:**
- [Installation](#installation) — CMake + build
- [Running Simulations](#running-simulations) — CLI commands
- [Configuration](#configuration) — Parameters and tuning
- [Output Files](#output-files) — Results format
- [Visualization](#visualization) — ParaView setup
- [Validation Cases](#validation-cases) — Test benchmarks
- [Troubleshooting](#troubleshooting) — Common issues

**See also:**
- Learning & theory → [COURSE.md](COURSE.md)
- Implementation details → [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md)
- Project overview → [README.md](README.md)

---

## Installation

### Prerequisites

- **CMake** ≥ 3.16
- **C++17 compiler** (g++ 7+, clang 6+, MSVC 2020+)
- **Eigen3** ≥ 3.3 (included/auto-fetched)
- **Optional:** CUDA Toolkit 11.0+ (GPU acceleration, currently optional)

### Build Steps

```bash
# Navigate to repository
cd cfd-engine

# Create build directory
mkdir -p build
cd build

# Configure with CMake (Release optimization recommended)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build
make -j4 # Use 4 parallel jobs

# Verify (executable created)
ls -lh cfd
```

### Build Options

```bash
# GPU acceleration (if CUDA available)
cmake -DCFD_BACKEND=CUDA ..

# CPU only (default)
cmake -DCFD_BACKEND=CPP ..

# Debugging (symbols + no optimization)
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Performance profiling
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-g" ..
```

### Troubleshooting Build

**"Eigen3 not found":**
```bash
# Install Eigen3
sudo apt-get install libeigen3-dev # Linux
brew install eigen # macOS
```

**"CUDA not found" (if needed):**
```bash
# Set CUDA path
cmake -DCUD A_TOOLKIT_ROOT_DIR=/usr/local/cuda ..
```

---

## Running Simulations

### Basic Commands

```bash
# Default: Lid cavity (Re=100)
./cfd

# Show all options
./cfd --help

# Lid cavity at different Reynolds
./cfd --re 1000 --nx 128 --ny 128

# Animation (50 frames)
./cfd --animate --numTimeSteps 50 --name my_sim

# Energy transport test
./cfd --case energy_transport

# Turbulent cavity (Re=10000)
./cfd --case high_re_cavity --re 10000 --nx 64 --ny 64

# Channel flow
./cfd --case channel_flow --nx 128 --ny 64
```

### CLI Reference

| Flag | Default | Description |
|------|---------|-------------|
| `--case` | `lid_cavity` | Test case to run |
| `--re` | `100` | Reynolds number |
| `--nx` | `64` | Grid points (x) |
| `--ny` | `64` | Grid points (y) |
| `--animate` | - | Generate animation |
| `--numTimeSteps` | `100` | Animation frame count |
| `--name` | `output` | Output directory name |
| `--help` | - | Show this help |

### Available Cases

| Case | Purpose | Typical Setup |
|------|---------|---------------|
| `lid_cavity` | Benchmark cavity with moving top | `--re 100 --nx 64` |
| `high_re_cavity` | Turbulent cavity | `--re 10000 --nx 64` (k-ε) |
| `channel_flow` | Turbulent channel | `--nx 128 --ny 64` (k-ε) |
| `energy_transport` | Heat diffusion test | Manufactured solution (MMS) |

---

## Configuration

### SolverConfig

**File:** `config/SolverConfig.hpp`

**Common parameters:**

```cpp
struct SolverConfig {
 double nu = 0.01; // Kinematic viscosity (adjusts Re)
 double dt = 0.1; // Time step (pseudo-transient)
 int maxIterations = 500; // Max SIMPLE iterations
 double alpha_u = 0.7; // U-velocity underrelax
 double alpha_v = 0.7; // V-velocity underrelax
 double alpha_p = 0.3; // Pressure underrelax
 double tolerance = 1e-4; // Residual tolerance
};
```

**Reynolds number:**
$$Re = \frac{U \times L}{\nu}$$

For unit domain with $U=1$: $Re = 1 / \nu$, so:
- Re=100: ν = 0.01
- Re=1000: ν = 0.001
- Re=10000: ν = 0.0001

### TurbulenceConfig

**For k-ε with wall functions:**

```cpp
auto cfg = config::TurbulenceConfig::withWallFunctions();
solver::StaggeredSIMPLE::solveSteady(
 mesh, u, v, p, nu, maxIter,
 alpha_u, alpha_v, alpha_p,
 verbosity, &cfg // ← Pass turbulence config
);
```

**Wall function activation:** Automatic when `y+ > 11.25`

### Tuning Relaxation Factors

**For convergence issues:**

```cpp
// If diverging, reduce pressure underrelaxation
alpha_p = 0.2 // Smaller = more stable, slower

// If oscillating, reduce velocity
alpha_u = 0.5
alpha_v = 0.5

// Typical ranges:
// Laminar: alpha_p ∈ [0.2, 0.4], alpha_u ∈ [0.6, 0.8]
// Turbulent: alpha_p ∈ [0.15, 0.3], alpha_u ∈ [0.4, 0.6]
```

---

## Output Files

### VTK/VTS Format (ParaView)

**Files:** `output/lid_cavity_{NNNNNNN}.vts` (7-digit frame numbering for up to 9.9M steps)

**Contains:**
- Velocity magnitude (colored field)
- Pressure field
- Turbulent viscosity (if turbulence enabled)
- Turbulent kinetic energy (if k-ε enabled)

**Open in ParaView:** Right-click → Open with ParaView

**Frame Numbering:** Frames are numbered with 7 digits (e.g., `0000000`, `0000001`, `9999999`) to support long simulations up to 10 million steps.

### CSV Files

**Residuals:** `output/csv/residuals.csv`
- Columns: iteration, mass_residual, u_residual, v_residual, k_residual, eps_residual, y_plus

**Field data:** `output/csv/field_*.csv`
- Snapshot of field at specific iteration

**Centerline:** `output/csv/centerline_u.csv`, `centerline_v.csv`
- For cavity: u and v profiles along centerlines

### Directory Structure

```
output/
├── lid_cavity_20260306_120000/ # Simulation run
│ ├── residuals.csv
│ ├── lid_cavity_0000000.vts # Frame 0 (7 digits)
│ ├── lid_cavity_0000001.vts # Frame 1
│ └── ...
├── csv/
│ ├── residuals.csv
│ └── centerline_u.csv
└── logs/
 └── simulation.log
```

### IO Error Handling

The solver implements a robust warning-and-continue error policy for output operations:

- **File write failures** (permission denied, disk full, etc.) are logged as warnings (`[IO WARNING]`) to stderr
- **Simulation continues** even if output fails, ensuring no loss of computational progress
- **CSV writer** throws exceptions for error conditions to maintain validation test compatibility
- **Format errors** during restart (VTS read) are logged as warnings; simulation continues from initial conditions if the restart file is unreadable

This approach prioritizes computational stability over output completeness.

---

## Visualization

### ParaView Setup

1. **Open ParaView**
 ```bash
 paraview
 ```

2. **Load VTS file:** File → Open → select `*.vts`

3. **Display options:**
 - Solid Color → Velocity Magnitude (colors cells)
 - Representation → Surface (smooth render)

4. **Add streamlines:**
 - Filters → Common → Streamlines
 - Set seed source, adjust density

5. **Animation:**
 - View → Animation View
 - Select time indicator (if PVD file exists)
 - Play

### Expected Results (Visual)

**Lid-driven cavity (Re=100):**
- Primary vortex center: ~(0.5, 0.6)
- Secondary vortex (bottom-left): visible
- Streamlines show circulation pattern

**High-Re cavity (Re=10000):**
- Primary vortex stronger, shifted toward corner
- Multiple secondary/tertiary vortices
- Sharp gradients near walls (use fine grid)

**Channel flow:**
- Parabolic velocity profile develops
- Higher velocities in center
- Zero at walls (no-slip)

---

## Validation Cases

### Test 1: Lid-Driven Cavity (Laminar)

**Reference:** Ghia et al. (1982), JCP 48:387-411

**Setup:**
```bash
./cfd --case lid_cavity --re 100 --nx 64 --ny 64
```

**Expected results:**
- SIMPLE converges in ~300 iterations
- Mass residual < 1e-10
- Centerline velocity u ≈ 0.2-0.3 (interior)
- Velocity profile matches Ghia reference (L2 error < 0.5)

**Command to extract centerline:**
```bash
grep "centerline" output/csv/residuals.csv | tail -1
```

### Test 2: Energy Transport (Manufactured Solution)

**Test equation:** 
$$T(x,y,t) = \cos(2\pi x)\cos(2\pi y) e^{-4\pi^2\alpha t}$$

**Setup:**
```bash
./cfd --case energy_transport
```

**Expected results:**
- L2 error < 0.01
- Relative error < 1%
- Solution decays exponentially
- Grid convergence: error ~ O(h²)

**Verify convergence:**
```bash
# Run on 3 grids
./cfd --case energy_transport --nx 32 --ny 32
./cfd --case energy_transport --nx 64 --ny 64
./cfd --case energy_transport --nx 128 --ny 128

# L2 error should halve when grid doubles (2nd-order)
```

### Test 3: Turbulent High-Re Cavity (k-ε)

**Setup:**
```bash
./cfd --case high_re_cavity --re 1000 --nx 64 --ny 64
```

**Expected results:**
- Convergence: ~500 iterations
- y+ in log-layer (30 < y+ < 50)
- Turbulent viscosity νt >> ν (100× in core)
- k and ε converge monotonically

**Check wall treatment:**
```bash
grep "Wall y+ diagnostics" output/logs/simulation.log
# Should show: max = 40-50, mean = 20-30 (in log-layer)
```

### Test 4: Turbulent Channel Flow

**Physics:** Fully developed turbulent channel at Re_τ = 5600

**Setup:**
```bash
./cfd --case channel_flow --nx 128 --ny 64
```

**Expected results:**
- Convergence: ~1500 iterations
- Parabolic velocity profile
- Peak turbulent kinetic energy near walls
- y+ > 30 (log-layer)

---

## Troubleshooting

### Simulation Diverges

**Symptoms:** Residuals increase, NaN values appear

**Solutions:**
1. Reduce underrelaxation:
 ```cpp
 alpha_p = 0.15 // Pressure
 alpha_u = 0.4 // Velocity
 ```

2. Decrease timestep (if transient):
 ```cpp
 dt = 0.05 // From 0.1
 ```

3. Use upwind only (no QUICK/MUSCL):
 ```cpp
 scheme = ConvectionScheme::Upwind
 ```

### ParaView: "Could not read data"

**Solution:** Ensure file path is correct
```bash
# Check file exists
ls -lh output/*.vts

# use full path if needed
paraview /path/to/output/cavit_00000.vts
```

### Grid Independence Not Achieved

**Problem:** Error doesn't decrease with refinement

**Check:**
1. Boundary conditions correct?
2. Sufficient SIMPLE iterations? (increase maxIterations)
3. Tolerance too loose? (decrease to 1e-6)

### Wall y+ Out of Range

**Problem:** y+ < 11 (viscous sublayer) but using wall functions

**Solution:** Refine near-wall mesh
```cpp
dy = Ly / (Ny * 2); // Double grid points
```

Or switch to wall-resolved (y+ < 1) without wall functions

### IO Warnings During Simulation

**Message:** `[IO WARNING] Could not open file for writing`

**Causes & Solutions:**
1. **Disk full** → Free up space or use a different output directory
2. **Permission denied** → Check write permissions on output directory:
 ```bash
 chmod 755 output/
 ```
3. **Path too long** → Use shorter directory paths (max ~260 chars on Windows)
4. **Invalid path** → Ensure output directory exists:
 ```bash
 mkdir -p output/results
 ```

**Impact:** Warnings do not stop the simulation; it continues computing but may skip frame writes. Check stderr (`2>` redirection) to capture warnings.

---

## Performance Tips

### Optimization

- **Compiler optimization:** `-DCMAKE_BUILD_TYPE=Release`
- **Parallel build:** `make -j8` (adjust to CPU count)
- **Large grids:** Use GPU (`-DCFD_BACKEND=CUDA`)

### Faster Convergence

- Use upwind (smaller matrix bandwidth)
- Reduce grid size for testing
- Increase underrelaxation (carefully)

### Memory Usage

**For different grid sizes:**
- 32×32: ~1 MB
- 64×64: ~10 MB
- 128×128: ~40 MB
- 256×256: ~150 MB (limit for SparseLU)

---

## References

- **Ghia et al. (1982):** Lid-driven cavity benchmark (JCP 48:387-411)
- **ParaView documentation:** https://www.paraview.org/
- **Course slides:** `course/slides/` (PDF lectures)

---

## See Also

- **[COURSE.md](COURSE.md)** — Theory and concepts
- **[TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md)** — Implementation details
- **[README.md](README.md)** — Quick overview

---

## Validation Benchmarks & Reference Data

### Available Benchmark Problems

Test cases validated using built-in reference data:
```bash
./build/cfd --case lid_cavity --re 100 # Ghia et al. (1982)
./build/cfd --case channel_flow --re 5600 # DNS data
./build/cfd --case energy_transport # MMS
```

### Included Benchmarks

**Lid-Driven Cavity (Laminar & Turbulent):**
- Ghia et al. (1982) at Re = 100, 400, 1000
- Extended to turbulent Re with k-ε model
- Centerline velocity profiles

**Channel Flow (Turbulent):**
- DNS reference data at Re_τ = 180, 395, 590
- Wall function range (11 < y+ < 300)
- Statistical profiles (mean velocity, RMS)

**Thermal Transport (MMS):**
- Manufactured solution with analytical error norms
- Grid convergence verification (O(h²) spatial, O(dt) temporal)

### High-Reynolds Cavity (Turbulent k-ε)

**Setup:** Square cavity (1m × 1m), Re = 500-2000+, Grid: 64×64 or finer
- Top wall: u = 1.0 m/s, other walls: no-slip

**Expected:**
- k peaks near walls, lower in core
- ε high near walls, moderate in core
- Wall y+: 11.25 < y+ < 300 (log-layer)
- Typical y+ ≈ 35-50 for Re=1000 on 64×64

**Run:**
```bash
./cfd --case high_re_cavity --re 1000
```

**Pass Criteria:**
- [Complete] Residual < 1e-6
- [Complete] y+ in [11, 300] 
- [Complete] k, ε > 0
- [Complete] ν_t/ν ≈ 10-100

### Turbulent Channel Flow

**Fully developed flow driven by pressure gradient:**
```bash
./cfd --case channel_flow --nx 128 --ny 64
```

**Expected:**
- Converges ~1500 iterations
- Parabolic velocity profile
- Peak k near walls
- y+ > 30 (log-layer)

