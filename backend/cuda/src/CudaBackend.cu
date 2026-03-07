// CudaBackend.cu - Clean CUDA backend implementation

#ifdef CFD_ENABLE_CUDA

#include "../include/CudaBackend.hpp"
#include <thrust/transform.h>
#include <thrust/inner_product.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#define CUDA_CHECK(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) \
        throw std::runtime_error(std::string("CUDA: ") + cudaGetErrorString(err)); \
} while(0)

#define CUSPARSE_CHECK(call) do { \
    cusparseStatus_t err = (call); \
    if (err != CUSPARSE_STATUS_SUCCESS) \
        throw std::runtime_error(std::string("cuSPARSE error: ") + std::to_string((int)err)); \
} while(0)

#define CUBLAS_CHECK(call) do { \
    cublasStatus_t err = (call); \
    if (err != CUBLAS_STATUS_SUCCESS) \
        throw std::runtime_error(std::string("cuBLAS error: ") + std::to_string((int)err)); \
} while(0)

namespace cfd::backend::cuda {

// ─── Global handles ─────────────────────────────────────────
static cusparseHandle_t g_cusparse = nullptr;
static cublasHandle_t   g_cublas   = nullptr;

static void initHandles() {
    if (!g_cusparse) {
        CUSPARSE_CHECK(cusparseCreate(&g_cusparse));
        CUBLAS_CHECK(cublasCreate(&g_cublas));
    }
}

// ─── Vector constructors ────────────────────────────────────
Vector::Vector(int n) : host_data(n, Scalar(0)), device_data(n, Scalar(0)) {}
Vector::Vector(int n, Scalar val) : host_data(n, val), device_data(n, val) {}

// Explicit copy/move — defined here so Thrust ops go through nvcc
Vector::Vector(const Vector& o)
    : host_data(o.host_data), device_data(o.device_data),
      host_is_dirty(o.host_is_dirty) {}

Vector::Vector(Vector&& o)
    : host_data(std::move(o.host_data)), device_data(std::move(o.device_data)),
      host_is_dirty(o.host_is_dirty) {}

Vector& Vector::operator=(const Vector& o) {
    if (this != &o) {
        host_data    = o.host_data;
        device_data  = o.device_data;
        host_is_dirty = o.host_is_dirty;
    }
    return *this;
}

Vector& Vector::operator=(Vector&& o) {
    if (this != &o) {
        host_data    = std::move(o.host_data);
        device_data  = std::move(o.device_data);
        host_is_dirty = o.host_is_dirty;
    }
    return *this;
}

// Explicit SparseMatrix copy/move
SparseMatrix::SparseMatrix(const SparseMatrix& o)
    : nrows(o.nrows), ncols(o.ncols), nnz(o.nnz),
      row_offsets(o.row_offsets), col_indices(o.col_indices), values(o.values) {}

SparseMatrix::SparseMatrix(SparseMatrix&& o)
    : nrows(o.nrows), ncols(o.ncols), nnz(o.nnz),
      row_offsets(std::move(o.row_offsets)),
      col_indices(std::move(o.col_indices)),
      values(std::move(o.values)) {}

SparseMatrix& SparseMatrix::operator=(const SparseMatrix& o) {
    if (this != &o) {
        nrows = o.nrows; ncols = o.ncols; nnz = o.nnz;
        row_offsets = o.row_offsets;
        col_indices = o.col_indices;
        values      = o.values;
    }
    return *this;
}

SparseMatrix& SparseMatrix::operator=(SparseMatrix&& o) {
    if (this != &o) {
        nrows = o.nrows; ncols = o.ncols; nnz = o.nnz;
        row_offsets = std::move(o.row_offsets);
        col_indices = std::move(o.col_indices);
        values      = std::move(o.values);
    }
    return *this;
}

Vector Vector::Zero(int n) { return Vector(n, Scalar(0)); }
Vector Vector::Ones(int n) { return Vector(n, Scalar(1)); }

void Vector::resize(int n) {
    host_data.resize(n, Scalar(0));
    device_data.resize(n, Scalar(0));
}

void Vector::uploadToDevice() const {
    device_data.assign(host_data.begin(), host_data.end());
}

void Vector::downloadToHost() const {
    // Use cudaMemcpy to copy device→host in a const-safe way
    size_t n = device_data.size();
    if (n == 0) return;
    const_cast<thrust::host_vector<Scalar>&>(host_data).assign(n, Scalar(0));
    CUDA_CHECK(cudaMemcpy(
        const_cast<thrust::host_vector<Scalar>&>(host_data).data(),
        thrust::raw_pointer_cast(device_data.data()),
        n * sizeof(Scalar),
        cudaMemcpyDeviceToHost));
    const_cast<bool&>(host_is_dirty) = false;
}

Scalar* Vector::devPtr() {
    uploadToDevice();
    return thrust::raw_pointer_cast(device_data.data());
}

const Scalar* Vector::devPtr() const {
    uploadToDevice();
    return thrust::raw_pointer_cast(device_data.data());
}

// ─── Vector arithmetic (host-side for simplicity & correctness) ─
Vector Vector::operator+(const Vector& o) const {
    syncToHost(); o.syncToHost();
    Vector result(size());
    for (int i = 0; i < size(); ++i) result.host_data[i] = host_data[i] + o.host_data[i];
    return result;
}

Vector Vector::operator-(const Vector& o) const {
    syncToHost(); o.syncToHost();
    Vector result(size());
    for (int i = 0; i < size(); ++i) result.host_data[i] = host_data[i] - o.host_data[i];
    return result;
}

Vector Vector::operator*(Scalar s) const {
    syncToHost();
    Vector result(size());
    for (int i = 0; i < size(); ++i) result.host_data[i] = host_data[i] * s;
    return result;
}

Vector& Vector::operator+=(const Vector& o) {
    syncToHost(); o.syncToHost();
    for (int i = 0; i < size(); ++i) host_data[i] += o.host_data[i];
    return *this;
}

Vector& Vector::operator-=(const Vector& o) {
    syncToHost(); o.syncToHost();
    for (int i = 0; i < size(); ++i) host_data[i] -= o.host_data[i];
    return *this;
}

Vector& Vector::operator*=(Scalar s) {
    syncToHost();
    for (int i = 0; i < size(); ++i) host_data[i] *= s;
    return *this;
}

// ─── Linear algebra ─────────────────────────────────────────
Scalar Vector::dot(const Vector& o) const {
    syncToHost(); o.syncToHost();
    Scalar sum = 0;
    for (int i = 0; i < size(); ++i) sum += host_data[i] * o.host_data[i];
    return sum;
}

Scalar Vector::squaredNorm() const { return dot(*this); }
Scalar Vector::norm()        const { return std::sqrt(squaredNorm()); }

void Vector::setZero() {
    std::fill(host_data.begin(), host_data.end(), Scalar(0));
}

void Vector::setOnes() {
    std::fill(host_data.begin(), host_data.end(), Scalar(1));
}

// ─── SparseMatrix::setFromTriplets (CPU sort → upload) ──────
void SparseMatrix::setFromTriplets(const Triplet* triplets, int count) {
    nnz = count;
    if (count == 0) {
        std::vector<int> h_offsets(nrows + 1, 0);
        row_offsets = h_offsets;
        col_indices.clear();
        values.clear();
        return;
    }

    // Sort triplets by row then col (on host)
    std::vector<Triplet> sorted(triplets, triplets + count);
    std::sort(sorted.begin(), sorted.end(), [](const Triplet& a, const Triplet& b) {
        return a.row < b.row || (a.row == b.row && a.col < b.col);
    });

    // Build CSR on host
    std::vector<int>    h_row_offsets(nrows + 1, 0);
    std::vector<int>    h_col_indices(count);
    std::vector<Scalar> h_values(count);

    for (int k = 0; k < count; ++k) {
        h_row_offsets[sorted[k].row + 1]++;
        h_col_indices[k] = sorted[k].col;
        h_values[k]      = sorted[k].value;
    }
    // Prefix sum
    for (int i = 0; i < nrows; ++i)
        h_row_offsets[i + 1] += h_row_offsets[i];

    // Upload to device
    row_offsets = h_row_offsets;
    col_indices = h_col_indices;
    values      = h_values;
}

// ─── SpMV via cuSPARSE ──────────────────────────────────────
Vector SparseMatrix::operator*(const Vector& x) const {
    if (nnz == 0) return Vector(nrows, Scalar(0));

    initHandles();

    // Ensure x is on device
    x.uploadToDevice();

    Vector y(nrows, Scalar(0));
    y.uploadToDevice();

    cusparseSpMatDescr_t matA;
    CUSPARSE_CHECK(cusparseCreateCsr(
        &matA, nrows, ncols, nnz,
        const_cast<int*>(thrust::raw_pointer_cast(row_offsets.data())),
        const_cast<int*>(thrust::raw_pointer_cast(col_indices.data())),
        const_cast<Scalar*>(thrust::raw_pointer_cast(values.data())),
        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_BASE_ZERO, CUDA_R_TYPE));

    cusparseDnVecDescr_t vecX, vecY;
    CUSPARSE_CHECK(cusparseCreateDnVec(
        &vecX, ncols, const_cast<Scalar*>(x.devPtr()), CUDA_R_TYPE));
    CUSPARSE_CHECK(cusparseCreateDnVec(
        &vecY, nrows, y.devPtr(), CUDA_R_TYPE));

    Scalar alpha = 1, beta = 0;
    size_t bufSize = 0;
    CUSPARSE_CHECK(cusparseSpMV_bufferSize(
        g_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha, matA, vecX, &beta, vecY, CUDA_R_TYPE,
        CUSPARSE_SPMV_ALG_DEFAULT, &bufSize));

    void* buf = nullptr;
    CUDA_CHECK(cudaMalloc(&buf, bufSize > 0 ? bufSize : 1));

    CUSPARSE_CHECK(cusparseSpMV(
        g_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha, matA, vecX, &beta, vecY, CUDA_R_TYPE,
        CUSPARSE_SPMV_ALG_DEFAULT, buf));

    CUDA_CHECK(cudaFree(buf));
    cusparseDestroySpMat(matA);
    cusparseDestroyDnVec(vecX);
    cusparseDestroyDnVec(vecY);

    // Bring result back to host
    y.host_is_dirty = true;
    y.syncToHost();

    return y;
}

// ─── Helper: SpMV wrapper ────────────────────────────────────
static Vector spmv(const SparseMatrix& A, const Vector& x) {
    return A * x;
}

// ─── CG Solver ───────────────────────────────────────────────
Vector solveCG(const SparseMatrix& A, const Vector& b, int maxIter, Scalar tol, int& iters) {
    int n = A.nrows;
    Vector x(n, Scalar(0));
    Vector r = b;
    Vector p = r;
    Scalar rsold = r.dot(r);

    for (iters = 0; iters < maxIter; ++iters) {
        Vector Ap = spmv(A, p);
        Scalar pAp = p.dot(Ap);
        if (std::abs(pAp) < Scalar(1e-16)) break;
        Scalar alpha = rsold / pAp;
        x += alpha * p;
        r -= alpha * Ap;
        Scalar rsnew = r.dot(r);
        if (std::sqrt(rsnew) < tol) { ++iters; break; }
        Scalar beta = rsnew / rsold;
        p = r + beta * p;
        rsold = rsnew;
    }
    return x;
}

// ─── BiCGSTAB Solver ─────────────────────────────────────────
Vector solveBiCGSTAB(const SparseMatrix& A, const Vector& b, int maxIter, Scalar tol, int& iters) {
    int n = A.nrows;
    Vector x(n, Scalar(0));
    Vector r = b;
    Vector r_hat = r;
    Scalar rho_prev = 1, alpha = 1, omega = 1;
    Vector v(n, Scalar(0)), p(n, Scalar(0));

    for (iters = 0; iters < maxIter; ++iters) {
        Scalar rho = r_hat.dot(r);
        if (std::abs(rho) < Scalar(1e-16)) break;
        if (iters == 0) {
            p = r;
        } else {
            Scalar beta = (rho / rho_prev) * (alpha / omega);
            p = r + beta * (p - omega * v);
        }
        v = spmv(A, p);
        Scalar r_hat_v = r_hat.dot(v);
        if (std::abs(r_hat_v) < Scalar(1e-16)) break;
        alpha = rho / r_hat_v;
        Vector s = r - alpha * v;
        if (s.norm() < tol) { x += alpha * p; ++iters; break; }
        Vector t = spmv(A, s);
        Scalar tt = t.dot(t);
        omega = (tt > Scalar(1e-16)) ? t.dot(s) / tt : Scalar(1);
        x += alpha * p + omega * s;
        r = s - omega * t;
        if (r.norm() < tol) { ++iters; break; }
        rho_prev = rho;
    }
    return x;
}

// ─── LU Solver (fallback to BiCGSTAB on GPU) ─────────────────
Vector solveLU(const SparseMatrix& A, const Vector& b) {
    int iters = 0;
    return solveBiCGSTAB(A, b, 10000, Scalar(1e-10), iters);
}

// ─── Arithmetic operators (defined here, NOT in header) ──────
// (friend declarations allow finding these via ADL)

}  // namespace cfd::backend::cuda

#endif  // CFD_ENABLE_CUDA
