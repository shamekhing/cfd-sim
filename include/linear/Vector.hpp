#pragma once

#include "../../backend/BackendDispatcher.hpp"

namespace linear {

// Scalar type
#ifndef CFD_SCALAR_FLOAT
#define CFD_SCALAR_FLOAT 0
#endif

#if CFD_SCALAR_FLOAT
using Scalar = float;
#else
using Scalar = double;
#endif

class Vector {
protected:
    ::BackendImpl::Vector data_;

public:
    Vector() = default;
    explicit Vector(int n);
    Vector(const Vector& other) = default;
    Vector& operator=(const Vector& other) = default;

    int size() const { return data_.size(); }
    
    Scalar& operator[](int i) { return data_(i); }
    Scalar operator[](int i) const { return data_(i); }
    
    Scalar operator()(int i) const { return data_(i); }
    Scalar& operator()(int i) { return data_(i); }
    
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(Scalar scalar) const;
    
    // Inline hot operations (frequently called in solvers/assembly)
    inline Scalar dot(const Vector& other) const {
        return data_.dot(other.data_);
    }
    
    inline Scalar norm() const {
        return data_.norm();
    }
    
    inline Scalar norm2() const {
        return data_.squaredNorm();
    }
    
    void setZero();
    void setOnes();
    
    // Support Scalar * Vector (left scalar multiplication)
    friend Vector operator*(Scalar scalar, const Vector& v) {
        return v * scalar;
    }
    
    // In-place arithmetic operations (inlined for performance)
    inline Vector& operator+=(const Vector& other) {
        data_ += other.data_;
        return *this;
    }
    
    inline Vector& operator-=(const Vector& other) {
        data_ -= other.data_;
        return *this;
    }
    
    inline Vector& operator*=(Scalar scalar) {
        data_ *= scalar;
        return *this;
    }
    
    inline ::BackendImpl::Vector& raw() { return data_; }
    inline const ::BackendImpl::Vector& raw() const { return data_; }
    
    inline void setRawData(const ::BackendImpl::Vector& v) { data_ = v; }
    inline void setRawData(::BackendImpl::Vector&& v) { data_ = std::move(v); }  // Move overload

private:
    // Private constructor for optimized operator returns (enables RVO/move elision)
    explicit Vector(const ::BackendImpl::Vector& data) : data_(data) {}
    explicit Vector(::BackendImpl::Vector&& data) : data_(std::move(data)) {}
};

// Optimized helper: direct backend initialization (avoids double setZero call)
inline Vector zeros(int n) {
    Vector v;
    v.raw() = ::BackendImpl::Vector::Zero(n);
    return v;
}

inline Vector ones(int n) {
    Vector v;
    v.raw() = ::BackendImpl::Vector::Ones(n);
    return v;
}

}  // namespace linear
