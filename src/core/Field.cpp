#include "../../include/core/Field.hpp"
#include "../../include/core/Indexing.hpp"

namespace core {

Field::Field(const Mesh& mesh)
{
    resize(mesh);
}

void Field::resize(const Mesh& mesh)
{
    Nx_ = mesh.Nx();
    Ny_ = mesh.Ny();
    data_ = linear::Vector(Nx_ * Ny_);
    data_.setZero();
}

linear::Scalar& Field::operator()(int i, int j)
{
    return data_(id(i, j, Nx_));
}

linear::Scalar Field::operator()(int i, int j) const
{
    return data_(id(i, j, Nx_));
}

linear::Scalar& Field::operator[](int index)
{
    return data_(index);
}

linear::Scalar Field::operator[](int index) const
{
    return data_(index);
}

void Field::setZero()
{
    data_.setZero();
}

}
