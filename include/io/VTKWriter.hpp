#pragma once

#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "ColorMap.hpp"

namespace io {

class VTKWriter {
public:

    static void write(
        const std::string& filename,
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v,
        const core::Field& p)
    {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "[IO WARNING] Could not open file for writing: " << filename << std::endl;
            return;
        }

        int Nx = mesh.Nx();
        int Ny = mesh.Ny();

        file << "# vtk DataFile Version 3.0\n";
        file << "CFD Output\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_GRID\n";
        file << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
        file << "POINTS " << Nx*Ny << " float\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                file << mesh.x(i) << " "
                     << mesh.y(j) << " 0\n";
            }
        }

        file << "\nPOINT_DATA " << Nx*Ny << "\n";

        // Velocity magnitude (most important for visualization!)
        file << "SCALARS velocity_magnitude float 1\n";
        file << "LOOKUP_TABLE default\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                double mag = std::sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j));
                file << mag << "\n";
            }
        }

        // Pressure scalar
        file << "\nSCALARS pressure float 1\n";
        file << "LOOKUP_TABLE default\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                file << p(i,j) << "\n";
            }
        }

        // Velocity vector
        file << "\nVECTORS velocity float\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                file << u(i,j) << " "
                     << v(i,j) << " 0\n";
            }
        }

        if (!file.good()) {
            std::cerr << "[IO WARNING] Error during write to file: " << filename << std::endl;
        }
        file.close();

        file.close();
    }

    /**
     * @brief Write fields including turbulence k, epsilon, and turbulent viscosity
     * 
     * @param filename Output VTK file path
     * @param mesh Computational mesh
     * @param u U-velocity field
     * @param v V-velocity field
     * @param p Pressure field
     * @param k Turbulent kinetic energy
     * @param epsilon Turbulent dissipation rate
     * @param nu_t Turbulent viscosity
     */
    static void writeWithTurbulence(
        const std::string& filename,
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v,
        const core::Field& p,
        const core::Field& k,
        const core::Field& epsilon,
        const core::Field& nu_t)
    {
        std::ofstream file(filename);

        int Nx = mesh.Nx();
        int Ny = mesh.Ny();

        file << "# vtk DataFile Version 3.0\n";
        file << "CFD Output with Turbulence\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_GRID\n";
        file << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
        file << "POINTS " << Nx*Ny << " float\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                file << mesh.x(i) << " "
                     << mesh.y(j) << " 0\n";
            }
        }

        file << "\nPOINT_DATA " << Nx*Ny << "\n";

        // Velocity magnitude
        file << "SCALARS velocity_magnitude float 1\n";
        file << "LOOKUP_TABLE default\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                double mag = std::sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j));
                file << mag << "\n";
            }
        }

        // Pressure scalar
        file << "\nSCALARS pressure float 1\n";
        file << "LOOKUP_TABLE default\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                file << p(i,j) << "\n";
            }
        }

        // Turbulent kinetic energy
        file << "\nSCALARS turbulent_k float 1\n";
        file << "LOOKUP_TABLE default\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                file << k(i,j) << "\n";
            }
        }

        // Turbulent dissipation rate
        file << "\nSCALARS turbulent_epsilon float 1\n";
        file << "LOOKUP_TABLE default\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                file << epsilon(i,j) << "\n";
            }
        }

        // Turbulent viscosity
        file << "\nSCALARS turbulent_viscosity float 1\n";
        file << "LOOKUP_TABLE default\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                file << nu_t(i,j) << "\n";
            }
        }

        // Velocity vector
        file << "\nVECTORS velocity float\n";

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                file << u(i,j) << " "
                     << v(i,j) << " 0\n";
            }
        }

        file.close();
    }

    /**
     * @brief Write with custom color scheme (Blue→Green→Red or other)
     * 
     * Generates scalar fields with embedded RGB lookup tables for enhanced visualization
     * 
     * @param filename Output VTK file path
     * @param mesh Computational mesh
     * @param u U-velocity field
     * @param v V-velocity field
     * @param p Pressure field
     * @param colorScheme Color scheme to apply
     */
    static void writeWithColorMap(
        const std::string& filename,
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v,
        const core::Field& p,
        ColorScheme colorScheme = ColorScheme::BLUE_GREEN_RED)
    {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "[IO WARNING] Could not open file for writing: " << filename << std::endl;
            return;
        }

        int Nx = mesh.Nx();
        int Ny = mesh.Ny();

        file << "# vtk DataFile Version 3.0\n";
        file << "CFD Output with Color Map\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_GRID\n";
        file << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
        file << "POINTS " << Nx*Ny << " float\n";

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                file << mesh.x(i) << " " << mesh.y(j) << " 0\n";
            }
        }

        file << "\nPOINT_DATA " << Nx*Ny << "\n";

        // Velocity magnitude with color map
        core::Field vel_mag(mesh);
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                vel_mag(i, j) = std::sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j));
            }
        }
        writeScalarWithColorMap(file, "velocity_magnitude", vel_mag, colorScheme);

        // Pressure with color map
        writeScalarWithColorMap(file, "pressure", p, colorScheme);

        // Velocity vector (not color-mapped, but included)
        file << "\nVECTORS velocity float\n";
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                file << u(i,j) << " " << v(i,j) << " 0\n";
            }
        }

        if (!file.good()) {
            std::cerr << "[IO WARNING] Error during write to file: " << filename << std::endl;
        }
        file.close();
    }

    /**
     * @brief Write vorticity and streamfunction fields with color scheme
     * 
     * @param filename Output VTK file path
     * @param mesh Computational mesh
     * @param vorticity Vorticity field (ω = ∂v/∂x - ∂u/∂y)
     * @param streamfunction Streamfunction field (ψ)
     * @param colorScheme Color scheme for visualization
     */
    static void writeVortexFields(
        const std::string& filename,
        const core::Mesh& mesh,
        const core::Field& vorticity,
        const core::Field& streamfunction,
        ColorScheme colorScheme = ColorScheme::BLUE_GREEN_RED)
    {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "[IO WARNING] Could not open vortex field file: " << filename << std::endl;
            return;
        }

        int Nx = mesh.Nx();
        int Ny = mesh.Ny();

        file << "# vtk DataFile Version 3.0\n";
        file << "Vortex Fields (Vorticity and Streamfunction)\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_GRID\n";
        file << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
        file << "POINTS " << Nx*Ny << " float\n";

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                file << mesh.x(i) << " " << mesh.y(j) << " 0\n";
            }
        }

        file << "\nPOINT_DATA " << Nx*Ny << "\n";

        // Vorticity magnitude with color map
        core::Field vort_mag(mesh);
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                vort_mag(i, j) = std::abs(vorticity(i, j));
            }
        }
        writeScalarWithColorMap(file, "vorticity_magnitude", vort_mag, colorScheme);

        // Streamfunction with color map
        writeScalarWithColorMap(file, "streamfunction", streamfunction, colorScheme);

        file.close();
    }

private:

    /**
     * @brief Helper: write a scalar field with embedded RGB color lookup table
     * 
     * @param file Output file stream
     * @param scalar_name Name of the scalar field
     * @param field Scalar field data
     * @param colorScheme Color scheme to apply
     */
    static void writeScalarWithColorMap(
        std::ofstream& file,
        const std::string& scalar_name,
        const core::Field& field,
        ColorScheme colorScheme)
    {
        // Generate lookup table
        std::vector<RGB> lut = ColorMap::generateLookupTable(colorScheme);
        
        // Get field range
        auto [min_val, max_val] = ColorMap::getFieldRange(field);
        
        // Map field values to color indices
        std::vector<unsigned char> indices = ColorMap::fieldToColorIndices(field, min_val, max_val);
        
        // Write scalar field header
        file << "SCALARS " << scalar_name << " unsigned_char 3\n";
        file << "LOOKUP_TABLE " << scalar_name << "_lut 256\n";
        
        // Write 256-entry lookup table
        for (int i = 0; i < 256; ++i) {
            unsigned char r = static_cast<unsigned char>(lut[i].r * 255.0f);
            unsigned char g = static_cast<unsigned char>(lut[i].g * 255.0f);
            unsigned char b = static_cast<unsigned char>(lut[i].b * 255.0f);
            file << static_cast<int>(r) << " "
                 << static_cast<int>(g) << " "
                 << static_cast<int>(b) << "\n";
        }
        
        // Write data: RGB triplets (one per point)
        file << "\n";
        for (size_t idx = 0; idx < indices.size(); ++idx) {
            unsigned char col_idx = indices[idx];
            const RGB& color = lut[col_idx];
            unsigned char r = static_cast<unsigned char>(color.r * 255.0f);
            unsigned char g = static_cast<unsigned char>(color.g * 255.0f);
            unsigned char b = static_cast<unsigned char>(color.b * 255.0f);
            file << static_cast<int>(r) << " "
                 << static_cast<int>(g) << " "
                 << static_cast<int>(b) << "\n";
        }
        file << "\n";
    }
};

}