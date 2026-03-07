#!/bin/bash

# CFD Engine - Run All Cases Script
# Builds and runs all cases with both Eigen and CUDA backends
# Usage: ./run_all_cases.sh

cd "$(dirname "$0")" || exit 1

echo "========================================"
echo "CFD Engine - Build and Run All Cases"
echo "========================================"
echo ""

# Build Eigen backend
echo "===> Building Eigen Backend..."
mkdir -p build_eigen
cd build_eigen
cmake -DCFD_BACKEND=EIGEN -DCMAKE_BUILD_TYPE=Release .. > /dev/null 2>&1
make -j4 > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ Eigen backend built successfully"
else
    echo "✗ Eigen backend build failed"
    exit 1
fi
cd ..
echo ""

# Build CUDA backend
echo "===> Building CUDA Backend..."
mkdir -p build_cuda
cd build_cuda
cmake -DCFD_BACKEND=CUDA -DCMAKE_BUILD_TYPE=Release .. > /dev/null 2>&1
make -j4 > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ CUDA backend built successfully"
else
    echo "✗ CUDA backend build failed"
    exit 1
fi
cd ..
echo ""

# Function to run all cases for a given backend
run_cases() {
    local backend=$1
    local build_dir=$2
    local output_dir=$3
    
    echo "========================================"
    echo "Running Cases with $backend Backend"
    echo "========================================"
    echo ""
    
    # Case 0: Energy Transport (MMS Validation)
    echo "[1/6] Energy Transport (MMS)..."
    ./$build_dir/cfd --case energy_transport --name energy_test --nx 64 --ny 64 --animate --numTimeSteps 100 --output $output_dir
    
    # Case 1: Lid-Driven Cavity Re=100
    echo "[2/6] Lid-Driven Cavity (Re=100)..."
    ./$build_dir/cfd --case lid_cavity --name lid100 --re 100 --nx 64 --ny 64 --animate --numTimeSteps 500 --output $output_dir
    
    # Case 2: Lid-Driven Cavity Re=1000
    echo "[3/6] Lid-Driven Cavity (Re=1000)..."
    ./$build_dir/cfd --case lid_cavity --name lid1000 --re 1000 --nx 128 --ny 128 --animate --numTimeSteps 1000 --output $output_dir
    
    # Case 3: High Reynolds Cavity (turbulent)
    echo "[4/6] High Reynolds Cavity (Re=10000)..."
    ./$build_dir/cfd --case high_re_cavity --name high_re --re 10000 --nx 64 --ny 64 --animate --numTimeSteps 500 --output $output_dir
    
    # Case 4: Turbulent Channel Flow
    echo "[5/6] Turbulent Channel Flow..."
    ./$build_dir/cfd --case channel_flow --name channel --nx 128 --ny 64 --animate --numTimeSteps 1500 --output $output_dir
    
    # Case 5: Cylinder Flow (Vortex Shedding)
    echo "[6/6] Cylinder Flow (Vortex Shedding)..."
    ./$build_dir/cfd --case cylinder_flow --name cylinder --re 100 --nx 96 --ny 96 --animate --numTimeSteps 2000 --output $output_dir
    
    echo ""
    echo "✓ All $backend cases completed"
    echo "  Results in $output_dir/"
    echo ""
}

run_main_case() {
    local backend=$1
    local build_dir=$2
    local output_dir=$3
    
    echo "========================================"
    echo "Running Main Case with $backend Backend"
    echo "========================================"
    echo ""

    # Case 0: Energy Transport (MMS Validation)
    echo "Flagship Cylinder Flow (Vortex Shedding)..."
        ./$build_dir/cfd --case cylinder_main --nx 200 --ny 100 --re 1000 --name flagship --animate --numTimeSteps 1500 --output $output_dir
    
    echo ""
    echo "✓ All $backend cases completed"
    echo "  Results in $output_dir/"
    echo ""
}

# Run cases with CUDA backend
#run_cases "CUDA" "build_cuda" "build_cuda/output"
run_main_case "CUDA" "build_cuda" "build_cuda/output"

# Run cases with Eigen backend
# run_cases "Eigen" "build_eigen" "build_eigen/output"

echo "========================================"
echo "All backends and cases completed!"
echo ""
echo "Results:"
echo "  Eigen: build_eigen/output/"
echo "  CUDA:  build_cuda/output/"
echo "========================================"