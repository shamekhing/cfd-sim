#pragma once

#include "../../backend/BackendDispatcher.hpp"
#include "Vector.hpp"
#include <vector>

namespace linear {

class SparseMatrix {
protected:
    ::BackendImpl::SparseMatrix data_;

public:
    SparseMatrix() = default;
    SparseMatrix(int rows, int cols);
    SparseMatrix(const SparseMatrix& other) = default;
    SparseMatrix& operator=(const SparseMatrix& other) = default;

    int rows() const { return data_.rows(); }
    int cols() const { return data_.cols(); }
    int nnz() const { return data_.nonZeros(); }
    
    void setFromTriplets(const ::BackendImpl::Triplet* triplets, int count);
    void makeCompressed();
    
    Vector operator*(const Vector& v) const;
    
    ::BackendImpl::SparseMatrix& raw() { return data_; }
    const ::BackendImpl::SparseMatrix& raw() const { return data_; }
};

class MatrixBuilder {
private:
    int rows_, cols_;
    std::vector<::BackendImpl::Triplet> triplets_;

public:
    MatrixBuilder(int rows, int cols) : rows_(rows), cols_(cols) {
        // Optimize reservation for typical CFD 2D stencil patterns
        // Standard 5-point Laplacian: ~5 neighbors per row
        int estimated_nnz = rows * 6;  // Slightly generous for boundary handling
        int max_nnz = rows * cols / 10;  // Cap at 10% sparsity for safety
        triplets_.reserve(std::min(estimated_nnz, max_nnz));
    }

    void add(int i, int j, Scalar value) {
        triplets_.emplace_back(i, j, value);
    }

    SparseMatrix build() {
        SparseMatrix A(rows_, cols_);
        A.setFromTriplets(triplets_.data(), triplets_.size());
        A.makeCompressed();
        return A;
    }

    void clear() {
        triplets_.clear();
    }
};

// Alias for backward compatibility with old assembly code
using SparseBuilder = MatrixBuilder;

}  // namespace linear
