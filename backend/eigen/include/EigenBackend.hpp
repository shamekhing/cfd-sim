#pragma once

// Eigen backend - ALL Eigen includes and implementations ONLY in this file
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <stdexcept>

// Enable parallel sparse operations
#ifndef EIGEN_NO_DEBUG
#define EIGEN_DONT_PARALLELIZE
#endif

namespace cfd::backend::eigen {

// Scalar type
#ifndef CFD_SCALAR_FLOAT
#define CFD_SCALAR_FLOAT 0
#endif

#if CFD_SCALAR_FLOAT
using Scalar = float;
#else
using Scalar = double;
#endif

// Base types
using DenseMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
using Triplet = Eigen::Triplet<Scalar>;

// ============================================================
// Vector: Wrapper around Eigen::VectorXd with compatible API
// ============================================================
struct Vector {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> data;

    Vector() = default;
    explicit Vector(int size) : data(size) {}
    Vector(const Vector&) = default;
    Vector(Vector&&) = default;
    Vector& operator=(const Vector&) = default;
    Vector& operator=(Vector&&) = default;

    // Static constructors
    static Vector Zero(int n) { Vector v; v.data = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(n); return v; }
    static Vector Ones(int n) { Vector v; v.data = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Ones(n); return v; }
    static Vector Constant(int n, Scalar val) { Vector v; v.data = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Constant(n, val); return v; }

    int size() const { return data.size(); }
    void resize(int n) { data.resize(n); }
    
    Scalar& operator[](int i) { return data[i]; }
    const Scalar& operator[](int i) const { return data[i]; }
    Scalar& operator()(int i) { return data(i); }
    const Scalar& operator()(int i) const { return data(i); }

    Vector operator+(const Vector& o) const { Vector r; r.data = data + o.data; return r; }
    Vector operator-(const Vector& o) const { Vector r; r.data = data - o.data; return r; }
    Vector operator*(Scalar s) const { Vector r; r.data = data * s; return r; }
    Vector operator/(Scalar s) const { Vector r; r.data = data / s; return r; }
    Vector operator-() const { Vector r; r.data = -data; return r; }

    Vector& operator+=(const Vector& o) { data += o.data; return *this; }
    Vector& operator-=(const Vector& o) { data -= o.data; return *this; }
    Vector& operator*=(Scalar s) { data *= s; return *this; }

    friend Vector operator*(Scalar s, const Vector& v) { return v * s; }

    Scalar dot(const Vector& o) const { return data.dot(o.data); }
    Scalar norm() const { return data.norm(); }
    Scalar squaredNorm() const { return data.squaredNorm(); }
    void setZero() { data.setZero(); }
    void setOnes() { data.setOnes(); }
};

// ============================================================
// SparseMatrix: Wrapper around Eigen::SparseMatrix<Scalar>
// ============================================================
struct SparseMatrix {
    Eigen::SparseMatrix<Scalar> data;

    SparseMatrix() = default;
    SparseMatrix(int r, int c) : data(r, c) {}
    SparseMatrix(const SparseMatrix&) = default;
    SparseMatrix(SparseMatrix&&) = default;
    SparseMatrix& operator=(const SparseMatrix&) = default;
    SparseMatrix& operator=(SparseMatrix&&) = default;

    int rows() const { return data.rows(); }
    int cols() const { return data.cols(); }
    int nonZeros() const { return data.nonZeros(); }

    // Compatible interface: accepts pointer+count, converts to iterators
    void setFromTriplets(const Triplet* triplets, int count) {
        data.setFromTriplets(triplets, triplets + count);
    }

    void makeCompressed() { data.makeCompressed(); }

    Vector operator*(const Vector& x) const {
        Vector result;
        result.data = data * x.data;
        return result;
    }
};

// ========== SOLVER IMPLEMENTATIONS ==========

// BiCGSTAB solver - iterative method for general matrices
inline Vector solveBiCGSTAB(const SparseMatrix& A, const Vector& b, int maxIter, Scalar tol, int& iters) {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<Scalar>> solver;
    solver.setMaxIterations(maxIter);
    solver.setTolerance(tol);
    solver.compute(A.data);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("BiCGSTAB: matrix decomposition failed");
    }
    
    auto x_eigen = solver.solve(b.data);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("BiCGSTAB: solve failed");
    }
    
    iters = solver.iterations();
    Vector x;
    x.data = x_eigen;
    return x;
}

// CG solver - for symmetric positive definite matrices
inline Vector solveCG(const SparseMatrix& A, const Vector& b, int maxIter, Scalar tol, int& iters) {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<Scalar>, Eigen::Lower|Eigen::Upper> solver;
    solver.setMaxIterations(maxIter);
    solver.setTolerance(tol);
    solver.compute(A.data);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("CG: matrix decomposition failed");
    }
    
    auto x_eigen = solver.solve(b.data);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("CG: solve failed");
    }
    
    iters = solver.iterations();
    Vector x;
    x.data = x_eigen;
    return x;
}

// LU solver - direct sparse LU decomposition
inline Vector solveLU(const SparseMatrix& A, const Vector& b) {
    Eigen::SparseLU<Eigen::SparseMatrix<Scalar>, Eigen::COLAMDOrdering<int>> solver;
    solver.compute(A.data);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("SparseLU: decomposition failed");
    }
    
    auto x_eigen = solver.solve(b.data);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("SparseLU: solve failed");
    }
    
    Vector x;
    x.data = x_eigen;
    return x;
}

}  // namespace cfd::backend::eigen
