#include "../../include/linear/Matrix.hpp"
#include "../../backend/BackendDispatcher.hpp"
#include <stdexcept>

namespace linear {

SparseMatrix::SparseMatrix(int rows, int cols) 
    : data_(::BackendImpl::SparseMatrix(rows, cols)) {}

void SparseMatrix::setFromTriplets(const ::BackendImpl::Triplet* triplets, int count) {
    data_.setFromTriplets(triplets, count);
}

void SparseMatrix::makeCompressed() {
    data_.makeCompressed();
}

Vector SparseMatrix::operator*(const Vector& v) const {
    Vector result;
    result.setRawData(data_ * v.raw());
    return result;
}

}  // namespace linear
