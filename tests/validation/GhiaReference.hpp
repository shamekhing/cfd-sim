#pragma once

#include <vector>
#include <cmath>
#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"

namespace validation {

// Ghia et al. (1982) JCP 48:387-411 - Lid-driven cavity Re=100
struct GhiaRe100 {
    static const std::vector<double> y_ref;
    static const std::vector<double> u_ref;
    static const std::vector<double> x_ref;
    static const std::vector<double> v_ref;
};

// Ghia et al. (1982) JCP 48:387-411 - Lid-driven cavity Re=1000
struct GhiaRe1000 {
    static const std::vector<double> y_ref;
    static const std::vector<double> u_ref;
    static const std::vector<double> x_ref;
    static const std::vector<double> v_ref;
};

// L2 error between computed centerline and Ghia reference (linear interpolation)
inline double centerlineErrorU(const core::Mesh& mesh, const core::Field& u,
                               const std::vector<double>& y_ref,
                               const std::vector<double>& u_ref)
{
    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();
    const int ix_mid = Nx / 2;
    double sum = 0.0;
    int count = 0;

    for (size_t k = 0; k < y_ref.size(); ++k)
    {
        double y = y_ref[k];
        double u_ghia = u_ref[k];

        for (int j = 0; j < Ny - 1; ++j)
        {
            if (mesh.y(j) <= y && y <= mesh.y(j + 1))
            {
                double y0 = mesh.y(j), y1 = mesh.y(j + 1);
                double u0 = u(ix_mid, j), u1 = u(ix_mid, j + 1);
                double denom = y1 - y0;
                if (denom < 1e-20) denom = 1e-20;
                double u_interp = u0 + (u1 - u0) * (y - y0) / denom;
                sum += (u_interp - u_ghia) * (u_interp - u_ghia);
                ++count;
                break;
            }
        }
    }
    return (count > 0) ? std::sqrt(sum / count) : 0.0;
}

inline double centerlineErrorV(const core::Mesh& mesh, const core::Field& v,
                               const std::vector<double>& x_ref,
                               const std::vector<double>& v_ref)
{
    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();
    const int jy_mid = Ny / 2;
    double sum = 0.0;
    int count = 0;

    for (size_t k = 0; k < x_ref.size(); ++k)
    {
        double x = x_ref[k];
        double v_ghia = v_ref[k];

        for (int i = 0; i < Nx - 1; ++i)
        {
            if (mesh.x(i) <= x && x <= mesh.x(i + 1))
            {
                double x0 = mesh.x(i), x1 = mesh.x(i + 1);
                double v0 = v(i, jy_mid), v1 = v(i + 1, jy_mid);
                double denom = x1 - x0;
                if (denom < 1e-20) denom = 1e-20;
                double v_interp = v0 + (v1 - v0) * (x - x0) / denom;
                sum += (v_interp - v_ghia) * (v_interp - v_ghia);
                ++count;
                break;
            }
        }
    }
    return (count > 0) ? std::sqrt(sum / count) : 0.0;
}

}
