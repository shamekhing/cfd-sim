// Scalar and vector field container
#pragma once

#include "../linear/Vector.hpp"
#include "Indexing.hpp"
#include "Mesh.hpp"

namespace core {

class Field {
private:
    int Nx_;
    int Ny_;
    linear::Vector data_;

public:
    Field() = default;

    explicit Field(const Mesh& mesh);

    void resize(const Mesh& mesh);

    int size() const { return data_.size(); }

    linear::Scalar& operator()(int i, int j);
    linear::Scalar operator()(int i, int j) const;

    linear::Scalar& operator[](int index);
    linear::Scalar operator[](int index) const;

    void setZero();

    linear::Vector& raw() { return data_; }
    const linear::Vector& raw() const { return data_; }
};

}
