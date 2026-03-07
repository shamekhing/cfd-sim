#pragma once

#include <utility>

namespace core {

// Cell-centered index (pressure, cell count Nx*Ny)
inline int id(int i, int j, int Nx)
{
    return i + j * Nx;
}

// U at vertical face: face (i,j) = east face of cell (i,j); i in [0..Nx], j in [0..Ny-1]
inline int idxU(int i, int j, int Nx)
{
    return i + j * (Nx + 1);
}

// V at horizontal face: face (i,j) = north face of cell (i,j); i in [0..Nx-1], j in [0..Ny]
inline int idxV(int i, int j, int Nx)
{
    return i + j * Nx;
}

inline std::pair<int,int> ij(int index, int Nx)
{
    int j = index / Nx;
    int i = index % Nx;
    return {i, j};
}

inline bool isInterior(int i, int j, int Nx, int Ny)
{
    return (i > 0 && i < Nx-1 && j > 0 && j < Ny-1);
}

inline bool isBoundary(int i, int j, int Nx, int Ny)
{
    return !isInterior(i, j, Nx, Ny);
}

// U-face (i,j): i in [0..Nx], j in [0..Ny-1]. Boundary: left i==0, right i==Nx, bottom j==0, top j==Ny-1 (lid).
inline bool isUBoundary(int i, int j, int Nx, int Ny)
{
    return i == 0 || i == Nx || j == 0 || j == Ny - 1;
}

// V-face (i,j): i in [0..Nx-1], j in [0..Ny]. Boundary: left i==0, right i==Nx-1, bottom j==0, top j==Ny.
inline bool isVBoundary(int i, int j, int Nx, int Ny)
{
    return i == 0 || i == Nx - 1 || j == 0 || j == Ny;
}

}
