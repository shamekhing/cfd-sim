#pragma once

#include <fstream>
#include <string>
#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"

namespace io {

class BenchmarkOutput {
public:
    // Write u(y) along vertical centerline (x = 0.5)
    // and v(x) along horizontal centerline (y = 0.5)
    static void writeCenterlines(
        const std::string& filename,
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();

        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "[IO WARNING] Could not open file for writing: " << filename << std::endl;
            return;
        }

        file << "x,y,u,v\n";

        // Vertical centerline: x fixed at center
        const int ix_mid = Nx / 2;
        for (int j = 0; j < Ny; ++j)
        {
            double x = mesh.x(ix_mid);
            double y = mesh.y(j);
            file << x << "," << y << "," << u(ix_mid, j) << "," << v(ix_mid, j) << "\n";
        }

        // Horizontal centerline: y fixed at center
        const int jy_mid = Ny / 2;
        for (int i = 0; i < Nx; ++i)
        {
            double x = mesh.x(i);
            double y = mesh.y(jy_mid);
            file << x << "," << y << "," << u(i, jy_mid) << "," << v(i, jy_mid) << "\n";
        }

        if (!file.good()) {
            std::cerr << "[IO WARNING] Error during write to file: " << filename << std::endl;
        }
        file.close();
    }

    // Write u(y) along vertical centerline only (for Ghia comparison: x=0.5)
    static void writeVerticalCenterlineU(
        const std::string& filename,
        const core::Mesh& mesh,
        const core::Field& u)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int ix_mid = Nx / 2;

        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "[IO WARNING] Could not open file for writing: " << filename << std::endl;
            return;
        }
        
        file << "y,u\n";

        for (int j = 0; j < Ny; ++j)
        {
            file << mesh.y(j) << "," << u(ix_mid, j) << "\n";
        }
        
        if (!file.good()) {
            std::cerr << "[IO WARNING] Error during write to file: " << filename << std::endl;
        }
        file.close();
    }

    // Write v(x) along horizontal centerline only (for Ghia comparison: y=0.5)
    static void writeHorizontalCenterlineV(
        const std::string& filename,
        const core::Mesh& mesh,
        const core::Field& v)
    {
        const int Nx = mesh.Nx();
        const int Ny = mesh.Ny();
        const int jy_mid = Ny / 2;

        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "[IO WARNING] Could not open file for writing: " << filename << std::endl;
            return;
        }
        
        file << "x,v\n";

        for (int i = 0; i < Nx; ++i)
        {
            file << mesh.x(i) << "," << v(i, jy_mid) << "\n";
        }
        
        if (!file.good()) {
            std::cerr << "[IO WARNING] Error during write to file: " << filename << std::endl;
        }
        file.close();
    }
};

}
