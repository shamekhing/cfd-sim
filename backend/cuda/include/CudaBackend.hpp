#pragma once

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <stdexcept>
#include <vector>

namespace cfd::backend::cuda {

// Scalar type
#ifndef CFD_SCALAR_FLOAT
#define CFD_SCALAR_FLOAT 0
#endif

#if CFD_SCALAR_FLOAT
using Scalar = float;
#define CUDA_R_TYPE CUDA_R_32F
#else
using Scalar = double;
#define CUDA_R_TYPE CUDA_R_64F
#endif

// Triplet for COO sparse matrix input
struct Triplet {
    int row, col;
    Scalar value;
    Triplet() = default;
    Triplet(int r, int c, Scalar v) : row(r), col(c), value(v) {}
};

// ============================================================
// Vector: host-mirror + device storage
// ALL copy/move operations defined in .cu to avoid Thrust in g++
// ============================================================
struct Vector {
    thrust::host_vector<Scalar>           host_data;
    mutable thrust::device_vector<Scalar> device_data;
    mutable bool host_is_dirty = false;

    Vector() = default;
    explicit Vector(int n);
    Vector(int n, Scalar val);

    // Explicitly declare copy/move so g++ doesn't instantiate Thrust internally
    Vector(const Vector& o);
    Vector(Vector&& o);
    Vector& operator=(const Vector& o);
    Vector& operator=(Vector&& o);

    static Vector Zero(int n);
    static Vector Ones(int n);

    int  size()   const { return (int)host_data.size(); }
    void resize(int n);

    // Host-side element access (wrapper uses these)
    Scalar& operator[](int i)       { return host_data[i]; }
    Scalar  operator[](int i) const { syncToHost(); return host_data[i]; }
    Scalar& operator()(int i)       { return host_data[i]; }
    Scalar  operator()(int i) const { syncToHost(); return host_data[i]; }

    // Arithmetic (all defined in .cu)
    Vector  operator+(const Vector& o) const;
    Vector  operator-(const Vector& o) const;
    Vector  operator*(Scalar s)        const;
    Vector& operator+=(const Vector& o);
    Vector& operator-=(const Vector& o);
    Vector& operator*=(Scalar s);

    friend Vector operator*(Scalar s, const Vector& v) { return v * s; }

    // Linear algebra (all defined in .cu)
    Scalar dot(const Vector& o) const;
    Scalar norm()        const;
    Scalar squaredNorm() const;
    void setZero();
    void setOnes();

    // Sync helpers (defined in .cu)
    void uploadToDevice() const;
    void downloadToHost() const;
    void syncToHost() const { if (host_is_dirty) downloadToHost(); }

    // Raw device pointer (uploads if needed, defined in .cu)
    Scalar*       devPtr();
    const Scalar* devPtr() const;
};

// ============================================================
// SparseMatrix: CSR format on GPU
// ALL copy/move operations defined in .cu to avoid Thrust in g++
// ============================================================
struct SparseMatrix {
    int nrows = 0, ncols = 0, nnz = 0;
    thrust::device_vector<int>    row_offsets;
    thrust::device_vector<int>    col_indices;
    thrust::device_vector<Scalar> values;

    SparseMatrix() = default;
    SparseMatrix(int r, int c) : nrows(r), ncols(c) {}

    // Explicit copy/move to ensure Thrust ops are compiled by nvcc
    SparseMatrix(const SparseMatrix& o);
    SparseMatrix(SparseMatrix&& o);
    SparseMatrix& operator=(const SparseMatrix& o);
    SparseMatrix& operator=(SparseMatrix&& o);

    int rows()     const { return nrows; }
    int cols()     const { return ncols; }
    int nonZeros() const { return nnz;   }

    void setFromTriplets(const Triplet* triplets, int count);
    void makeCompressed() {}

    Vector operator*(const Vector& x) const;
};

// ============================================================
// Solvers
// ============================================================
Vector solveBiCGSTAB(const SparseMatrix& A, const Vector& b, int maxIter, Scalar tol, int& iters);
Vector solveCG      (const SparseMatrix& A, const Vector& b, int maxIter, Scalar tol, int& iters);
Vector solveLU      (const SparseMatrix& A, const Vector& b);

}  // namespace cfd::backend::cuda
