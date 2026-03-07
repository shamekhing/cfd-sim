#pragma once

#include <iostream>
#include <filesystem>
#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "../include/solver/SIMPLE.hpp"
#include "../include/io/VTKWriter.hpp"

namespace cases {

inline void runSIMPLETest(
    const std::string& outputDir = "build/output",
    const std::string& name = "simple_test")
{
    const int Nx = 50;
    const int Ny = 50;

    core::Mesh mesh(Nx, Ny, 1.0, 1.0);

    core::Field u(mesh);
    core::Field v(mesh);
    core::Field p(mesh);

    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            u(i,j) = 0.0;
            v(i,j) = 0.0;
            p(i,j) = 0.0;
        }
    }

    const std::string vtkFile = outputDir + "/" + name + "/" + name + ".vtk";
    std::filesystem::create_directories(outputDir + "/" + name);
    solver::SIMPLE::solve(mesh, u, v, p, 0.01, 20);
    io::VTKWriter::write(vtkFile,mesh,u,v,p);
    std::cout << "SIMPLE solver finished. Output: " << vtkFile << "\n";
}

}