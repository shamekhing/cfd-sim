#pragma once

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"

namespace io {

class VTSReader {
public:
    // Read fields from VTS file and populate u, v, p
    static bool read(
        const std::string& filename,
        const core::Mesh& mesh,
        core::Field& u,
        core::Field& v,
        core::Field& p)
    {
        std::ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }

        int Nx = mesh.Nx();
        int Ny = mesh.Ny();
        int totalPoints = Nx * Ny;

        // Read entire file into memory for parsing
        std::stringstream buffer;
        buffer << file.rdbuf();
        std::string content = buffer.str();
        file.close();

        // Extract velocity vector data (contains u and v components)
        std::vector<double> u_data, v_data, p_data;
        
        if (!extractVelocityVectors(content, totalPoints, u_data, v_data)) {
            return false;
        }
        
        if (!extractPressure(content, totalPoints, p_data)) {
            return false;
        }

        // Populate Field objects (cell-centered data)
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int idx = j * Nx + i;
                u(i, j) = u_data[idx];
                v(i, j) = v_data[idx];
                p(i, j) = p_data[idx];
            }
        }

        return true;
    }

    // Read a named scalar field from a VTS file written by writeVTSFileScalar
    static bool readScalar(
        const std::string& filename,
        const core::Mesh& mesh,
        const std::string& scalarName,
        core::Field& scalar)
    {
        std::ifstream file(filename);
        if (!file.is_open()) return false;

        int Nx = mesh.Nx();
        int Ny = mesh.Ny();
        int totalPoints = Nx * Ny;

        std::stringstream buffer;
        buffer << file.rdbuf();
        std::string content = buffer.str();
        file.close();

        std::vector<double> s_data;
        if (!extractNamedScalar(content, scalarName, totalPoints, s_data)) return false;

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                scalar(i, j) = s_data[j * Nx + i];

        return true;
    }

private:
    // Extract velocity vector data (u, v components) using efficient streaming parser
    static bool extractVelocityVectors(const std::string& content, int totalPoints,
                                       std::vector<double>& u_data, std::vector<double>& v_data)
    {
        // Find velocity vector section without regex for better performance
        size_t velocityTag = content.find("Name=\"velocity\"");
        if (velocityTag == std::string::npos) {
            return false;
        }
        
        // Find the opening DataArray tag
        size_t dataStart = content.find(">", velocityTag);
        if (dataStart == std::string::npos) {
            return false;
        }
        dataStart++; // Move past the >
        
        // Find the closing DataArray tag
        size_t dataEnd = content.find("</DataArray>", dataStart);
        if (dataEnd == std::string::npos) {
            return false;
        }
        
        std::string velocityData = content.substr(dataStart, dataEnd - dataStart);
        
        // Parse velocity vectors (format: "u v 0\n" for each point)
        std::istringstream iss(velocityData);
        u_data.clear();
        v_data.clear();
        u_data.reserve(totalPoints);
        v_data.reserve(totalPoints);
        
        double u_val, v_val, w_val;
        int count = 0;
        while (iss >> u_val >> v_val >> w_val && count < totalPoints) {
            u_data.push_back(u_val);
            v_data.push_back(v_val);
            count++;
        }

        return count == totalPoints;
    }

    // Extract any named scalar field by name using efficient streaming parser
    static bool extractNamedScalar(const std::string& content, const std::string& fieldName,
                                   int totalPoints, std::vector<double>& data)
    {
        // Find the field by name without regex for better performance
        std::string nameTag = "Name=\"" + fieldName + "\"";
        size_t nameLoc = content.find(nameTag);
        if (nameLoc == std::string::npos) {
            return false;
        }
        
        // Find the opening DataArray tag
        size_t dataStart = content.find(">", nameLoc);
        if (dataStart == std::string::npos) {
            return false;
        }
        dataStart++; // Move past the >
        
        // Find the closing DataArray tag
        size_t dataEnd = content.find("</DataArray>", dataStart);
        if (dataEnd == std::string::npos) {
            return false;
        }
        
        std::string rawData = content.substr(dataStart, dataEnd - dataStart);
        std::istringstream iss(rawData);
        data.clear();
        data.reserve(totalPoints);
        double val;
        int count = 0;
        while (iss >> val && count < totalPoints) {
            data.push_back(val);
            count++;
        }
        return count == totalPoints;
    }

    // Extract pressure scalar data using efficient streaming parser
    static bool extractPressure(const std::string& content, int totalPoints,
                               std::vector<double>& p_data)
    {
        // Find pressure section without regex for better performance
        size_t pressureTag = content.find("Name=\"pressure\"");
        if (pressureTag == std::string::npos) {
            return false;
        }

        // Find the opening DataArray tag
        size_t dataStart = content.find(">", pressureTag);
        if (dataStart == std::string::npos) {
            return false;
        }
        dataStart++; // Move past the >
        
        // Find the closing DataArray tag
        size_t dataEnd = content.find("</DataArray>", dataStart);
        if (dataEnd == std::string::npos) {
            return false;
        }

        std::string pressureData = content.substr(dataStart, dataEnd - dataStart);
        
        // Parse pressure values
        std::istringstream iss(pressureData);
        p_data.clear();
        p_data.reserve(totalPoints);
        
        double p_val;
        int count = 0;
        while (iss >> p_val && count < totalPoints) {
            p_data.push_back(p_val);
            count++;
        }

        return count == totalPoints;
    }
};

}
