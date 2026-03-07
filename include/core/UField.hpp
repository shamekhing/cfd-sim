#pragma once

#include "../linear/Vector.hpp"
#include "Mesh.hpp"
#include "Field.hpp"
#include "Indexing.hpp"

namespace core {

// u-velocity at vertical faces; size (Nx+1)*Ny
// face (i,j) = east face of cell (i,j); i in [0..Nx], j in [0..Ny-1]
class UField {
public:
    UField() = default;
    explicit UField(const Mesh& mesh) {
        Nx_ = mesh.Nx();
        Ny_ = mesh.Ny();
        data_ = linear::Vector(mesh.totalUFaces());
        data_.setZero();
    }

    void resize(const Mesh& mesh) {
        Nx_ = mesh.Nx();
        Ny_ = mesh.Ny();
        data_ = linear::Vector(mesh.totalUFaces());
        data_.setZero();
    }

    int size() const { return data_.size(); }
    int Nx() const { return Nx_; }
    int Ny() const { return Ny_; }

    linear::Scalar& operator()(int i, int j) { return data_(idxU(i, j, Nx_)); }
    linear::Scalar operator()(int i, int j) const { return data_(idxU(i, j, Nx_)); }
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
                linear::Scalar uL = (i > 0) ? (*this)(i, j) : static_cast<linear::Scalar>(0.0);
                linear::Scalar uR = (i < Nx) ? (*this)(i + 1, j) : static_cast<linear::Scalar>(0.0);
                out(i, j) = 0.5 * (uL + uR);
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
