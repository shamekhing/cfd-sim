#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "../include/linear/Vector.hpp"

namespace io {

/**
 * @brief Simple CSV writer for exporting field data to CSV format
 */
class CSVWriter {
public:
    /**
     * @brief Write vector values to CSV file
     * 
     * @param filename Path to output CSV file
     * @param values Eigen vector to write
     * @param precision Decimal precision for floating point values
     */
    static void writeFieldToCSV(const std::string& filename,
                                const linear::Vector& values,
                                int precision = 8)
    {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open CSV file: " + filename);
        }

        file << std::scientific << std::setprecision(precision);
        for (int i = 0; i < values.size(); ++i) {
            file << values[i];
            if (i < values.size() - 1) file << "\n";
        }
        file << "\n";
        
        if (!file.good()) {
            throw std::runtime_error("Error writing to CSV file: " + filename);
        }
        file.close();
    }

    /**
     * @brief Write std::vector values to CSV file
     * 
     * @param filename Path to output CSV file
     * @param values Vector of values to write
     * @param precision Decimal precision for floating point values
     */
    static void writeFieldToCSV(const std::string& filename,
                                const std::vector<double>& values,
                                int precision = 8)
    {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open CSV file: " + filename);
        }

        file << std::scientific << std::setprecision(precision);
        for (size_t i = 0; i < values.size(); ++i) {
            file << values[i];
            if (i < values.size() - 1) file << "\n";
        }
        file << "\n";
        
        if (!file.good()) {
            throw std::runtime_error("Error writing to CSV file: " + filename);
        }
        file.close();
    }

    /**
     * @brief Write convergence data to CSV with headers
     * 
     * @param filename Path to output CSV file
     * @param headers Column headers
     * @param data 2D data (each row is a vector of values)
     * @param precision Decimal precision
     */
    static void writeConvergenceData(const std::string& filename,
                                     const std::vector<std::string>& headers,
                                     const std::vector<std::vector<double>>& data,
                                     int precision = 8)
    {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open CSV file: " + filename);
        }

        // Write headers
        for (size_t i = 0; i < headers.size(); ++i) {
            file << headers[i];
            if (i < headers.size() - 1) file << ",";
        }
        file << "\n";

        // Write data rows
        file << std::scientific << std::setprecision(precision);
        for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << row[i];
                if (i < row.size() - 1) file << ",";
            }
            file << "\n";
        }
        
        if (!file.good()) {
            throw std::runtime_error("Error writing to CSV file: " + filename);
        }
        file.close();
    }

    /**
     * @brief Write centerline profile to CSV
     * 
     * @param filename Path to output CSV file
     * @param coordinates Grid coordinates (e.g., y/Y for vertical line)
     * @param values Field values at coordinates
     * @param precision Decimal precision
     */
    static void writeProfile(const std::string& filename,
                             const std::vector<double>& coordinates,
                             const std::vector<double>& values,
                             int precision = 8)
    {
        if (coordinates.size() != values.size()) {
            throw std::runtime_error("Coordinate and value arrays must have same size");
        }

        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open CSV file: " + filename);
        }

        file << "x,y\n";  // Simple generic header
        file << std::scientific << std::setprecision(precision);
        for (size_t i = 0; i < coordinates.size(); ++i) {
            file << coordinates[i] << "," << values[i] << "\n";
        }
        
        if (!file.good()) {
            throw std::runtime_error("Error writing to CSV file: " + filename);
        }
        file.close();
    }
};

}  // namespace io