#include "../../include/linear/Vector.hpp"
#include <stdexcept>

namespace linear {

Vector::Vector(int n) : data_(::BackendImpl::Vector::Zero(n)) {}

Vector Vector::operator+(const Vector& other) const {
    return Vector(data_ + other.data_);
}

Vector Vector::operator-(const Vector& other) const {
    return Vector(data_ - other.data_);
}

Vector Vector::operator*(Scalar scalar) const {
    return Vector(data_ * scalar);
}

void Vector::setZero() {
    data_.setZero();
}

void Vector::setOnes() {
    data_.setOnes();
}

}  // namespace linear
