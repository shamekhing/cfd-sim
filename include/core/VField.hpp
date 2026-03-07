#pragma once

#include "../linear/Vector.hpp"
#include "Mesh.hpp"
#include "Field.hpp"
#include "Indexing.hpp"

namespace core {

// v-velocity at horizontal faces; size Nx*(Ny+1)
// face (i,j) = north face of cell (i,j); i in [0..Nx-1], j in [0..Ny]
class VField {
public:
    VField() = default;
    explicit VField(const Mesh& mesh) {
        Nx_ = mesh.Nx();
        Ny_ = mesh.Ny();
        data_ = linear::Vector(mesh.totalVFaces());
        data_.setZero();
    }

    void resize(const Mesh& mesh) {
        Nx_ = mesh.Nx();
        Ny_ = mesh.Ny();
        data_ = linear::Vector(mesh.totalVFaces());
        data_.setZero();
    }

    int size() const { return data_.size(); }
    int Nx() const { return Nx_; }
    int Ny() const { return Ny_; }

    linear::Scalar& operator()(int i, int j) { return data_(idxV(i, j, Nx_)); }
    linear::Scalar operator()(int i, int j) const { return data_(idxV(i, j, Nx_)); }
    linear::Scalar& operator[](int idx) { return data_(idx); }
    linear::Scalar operator[](int idx) const { return data_(idx); }

    void setZero() { data_.setZero(); }
    linear::Vector& raw() { return data_; }
    const linear::Vector& raw() const { return data_; }

    Field cellCenterView(const Mesh& mesh) const {
        Field out(mesh);
        int Nx = mesh.Nx();
        int Ny = mesh.Ny();
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                linear::Scalar vB = (j > 0) ? (*this)(i, j) : static_cast<linear::Scalar>(0.0);
                linear::Scalar vT = (j < Ny) ? (*this)(i, j + 1) : static_cast<linear::Scalar>(0.0);
                out(i, j) = 0.5 * (vB + vT);
            }
        }
        return out;
    }

private:
    int Nx_ = 0;
    int Ny_ = 0;
    linear::Vector data_;
};

}
