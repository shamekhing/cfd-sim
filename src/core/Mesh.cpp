#include "../../include/core/Mesh.hpp"

namespace core {

Mesh::Mesh(int Nx, int Ny, linear::Scalar Lx, linear::Scalar Ly)
    : Nx_(Nx), Ny_(Ny), Lx_(Lx), Ly_(Ly)
{
    dx_ = Lx_ / (Nx_ - 1);
    dy_ = Ly_ / (Ny_ - 1);
}

linear::Scalar Mesh::x(int i) const
{
    return i * dx_;
}

linear::Scalar Mesh::y(int j) const
{
    return j * dy_;
}

bool Mesh::isBoundary(int i, int j) const
{
    return (i == 0 || i == Nx_ - 1 ||
            j == 0 || j == Ny_ - 1);
}

bool Mesh::isInterior(int i, int j) const
{
    return !isBoundary(i, j);
}

}
