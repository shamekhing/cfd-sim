#pragma once

#include "../linear/Vector.hpp"

namespace core {

class Mesh {
private:
    int Nx_;
    int Ny_;
    linear::Scalar Lx_;
    linear::Scalar Ly_;
    linear::Scalar dx_;
    linear::Scalar dy_;

public:
    Mesh(int Nx, int Ny, linear::Scalar Lx, linear::Scalar Ly);

    int Nx() const { return Nx_; }
    int Ny() const { return Ny_; }

    int totalCells() const { return Nx_ * Ny_; }
    int totalUFaces() const { return (Nx_ + 1) * Ny_; }
    int totalVFaces() const { return Nx_ * (Ny_ + 1); }

    linear::Scalar Lx() const { return Lx_; }
    linear::Scalar Ly() const { return Ly_; }

    linear::Scalar dx() const { return dx_; }
    linear::Scalar dy() const { return dy_; }

    linear::Scalar x(int i) const;
    linear::Scalar y(int j) const;

    bool isLeftBoundary(int i) const { return i == 0; }
    bool isRightBoundary(int i) const { return i == Nx_ - 1; }
    bool isBottomBoundary(int j) const { return j == 0; }
    bool isTopBoundary(int j) const { return j == Ny_ - 1; }

    bool isBoundary(int i, int j) const;
    bool isInterior(int i, int j) const;
};

}
