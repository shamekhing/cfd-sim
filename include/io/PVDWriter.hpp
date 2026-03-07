#pragma once

#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <cstdlib>
#include <future>
#include <mutex>
#include <optional>
#include <regex>
#include "../include/core/Mesh.hpp"
#include "../include/core/Field.hpp"
#include "ColorMap.hpp"

namespace io {

class PVDWriter {
private:
    struct TimestepInfo {
        int step;
        double time;
        std::string filename;
    };
    
    std::vector<TimestepInfo> timesteps;
    std::string outputDir;
    std::string baseName;
    std::string pvdOutputDir;
    int verbosity;
    std::string absolutePath;
    ColorScheme colorScheme;
    bool useColorMap;
    
    // Async I/O support
    std::optional<std::future<void>> pendingWrite;
    std::mutex timestepsMutex;
    
    // Resume support
    int startFrame;
    std::string lastFrameFilename;

public:
    PVDWriter(const std::string& dirName = "animation", const std::string& name = "output", int verb = 1,
              const std::string& pvdDir = ".", ColorScheme cs = ColorScheme::BLUE_GREEN_RED, bool useColors = true)
        : outputDir(dirName), baseName(name), pvdOutputDir(pvdDir), verbosity(verb), 
          colorScheme(cs), useColorMap(useColors), startFrame(0), lastFrameFilename("")
    {
        // Create output directory if it doesn't exist
        std::filesystem::create_directories(outputDir);
        std::filesystem::create_directories(pvdOutputDir);
        
        // Get absolute path for PVD file references
        absolutePath = std::filesystem::absolute(outputDir).string();
        
        // Detect existing frames for resume
        detectExistingFrames();
        
        // Load existing PVD entries if resuming
        if (startFrame > 0) {
            loadExistingPVD();
        }
    }
    
    // Get the starting frame number (0 for new, N for resume from frame N)
    int getStartFrame() const { return startFrame; }
    
    // Get the filename of the last frame (for loading field state)
    std::string getLastFrameFilename() const { return lastFrameFilename; }

    // Write a single VTS (XML VTK) file for a specific timestep (async)
    void writeTimestep(
        int step,
        double time,
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v,
        const core::Field& p)
    {
        // Wait for previous async write to complete (single background thread)
        if (pendingWrite.has_value()) {
            pendingWrite->get();
            pendingWrite.reset();
        }

        // Generate filename with zero-padded step number (.vts = XML StructuredGrid)
        std::ostringstream oss;
        oss << outputDir << "/" << baseName << "_" 
            << std::setfill('0') << std::setw(7) << step << ".vts";
        std::string filename = oss.str();

        // Create relative filename for PVD metadata (7-digit padding for long runs)
        std::string relativeVtkName = baseName + "_" + 
            std::string(7 - std::to_string(step).length(), '0') + std::to_string(step) + ".vts";

        // Make deep copies of field data for async thread (avoid race conditions)
        // Use move semantics in lambda capture to avoid extra copy during async launch
        core::Field u_copy = u;
        core::Field v_copy = v;
        core::Field p_copy = p;

        // Launch async file write in background thread with move-captured fields
        pendingWrite = std::async(std::launch::async, 
            [this, filename, step, time, relativeVtkName, 
             u_copy = std::move(u_copy), v_copy = std::move(v_copy), p_copy = std::move(p_copy), &mesh]() mutable {
                // Write VTS file in background
                if (useColorMap) {
                    writeVTSFileWithColors(filename, mesh, u_copy, v_copy, p_copy, colorScheme);
                } else {
                    writeVTSFile(filename, mesh, u_copy, v_copy, p_copy);
                }
                
                // Thread-safe append to timesteps vector
                {
                    std::lock_guard<std::mutex> lock(timestepsMutex);
                    timesteps.push_back({step, time, relativeVtkName});
                }

                // Optional verbose logging (only occasional)
                if (verbosity > 1 && (step % 10 == 0 || step == 0)) {
                    std::lock_guard<std::mutex> lock(timestepsMutex);
                    std::cout << "  Wrote: " << filename << std::endl;
                }
            });
    }

    // Overload: write a single primary scalar field with a custom name (e.g. "temperature")
    // velocity u/v are still written but scalar is set as the default display variable
    void writeTimestepScalar(
        int step,
        double time,
        const core::Mesh& mesh,
        const core::Field& scalar,
        const std::string& scalarName,
        const core::Field& u,
        const core::Field& v)
    {
        if (pendingWrite.has_value()) {
            pendingWrite->get();
            pendingWrite.reset();
        }

        std::ostringstream oss;
        oss << outputDir << "/" << baseName << "_"
            << std::setfill('0') << std::setw(7) << step << ".vts";
        std::string filename = oss.str();

        // 7-digit padding for long runs
        std::string relativeVtkName = baseName + "_" +
            std::string(7 - std::to_string(step).length(), '0') + std::to_string(step) + ".vts";

        core::Field s_copy = scalar;
        core::Field u_copy = u;
        core::Field v_copy = v;

        // Use move semantics in lambda capture to avoid extra copies
        pendingWrite = std::async(std::launch::async,
            [this, filename, step, time, relativeVtkName, scalarName,
             s_copy = std::move(s_copy), u_copy = std::move(u_copy), v_copy = std::move(v_copy), &mesh]() mutable {
                writeVTSFileScalar(filename, mesh, s_copy, scalarName, u_copy, v_copy);
                {
                    std::lock_guard<std::mutex> lock(timestepsMutex);
                    timesteps.push_back({step, time, relativeVtkName});
                }
                if (verbosity > 1 && (step % 10 == 0 || step == 0)) {
                    std::lock_guard<std::mutex> lock(timestepsMutex);
                    std::cout << "  Wrote: " << filename << std::endl;
                }
            });
    }

    // Wait for all pending async writes to complete
    void finalize()
    {
        if (pendingWrite.has_value()) {
            pendingWrite->get();
            pendingWrite.reset();
        }
    }

    // Finalize and write the PVD metafile
    void writePVD()
    {
        std::string pvdFilename = (std::filesystem::path(pvdOutputDir) / (baseName + ".pvd")).string();

        std::ofstream pvdFile(pvdFilename);
        if (!pvdFile.is_open()) {
            std::cerr << "[IO WARNING] Could not open PVD file for writing: " << pvdFilename << std::endl;
            return;
        }
        
        pvdFile << "<?xml version=\"1.0\"?>\n";
        pvdFile << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
        pvdFile << "  <Collection>\n";

        for (const auto& ts : timesteps) {
            std::filesystem::path vtsPath = std::filesystem::path(outputDir) / ts.filename;
            std::error_code ec;
            std::filesystem::path relPath = std::filesystem::relative(vtsPath, std::filesystem::path(pvdOutputDir), ec);
            std::string fullPath = ec ? vtsPath.generic_string() : ("./" + relPath.generic_string());
            pvdFile << "    <DataSet timestep=\"" << std::fixed << std::setprecision(6) 
                    << ts.time << "\" file=\"" << fullPath << "\"/>\n";
        }

        pvdFile << "  </Collection>\n";
        pvdFile << "</VTKFile>\n";

        if (!pvdFile.good()) {
            std::cerr << "[IO WARNING] Error during write to PVD file: " << pvdFilename << std::endl;
        }
        pvdFile.close();

        if (verbosity > 0) {
            std::cout << "  Wrote PVD collection: " << pvdFilename << std::endl;
            std::cout << "  Total timesteps: " << timesteps.size() << std::endl;
        }
    }

    // Clear timestep history for new simulation
    void reset()
    {
        timesteps.clear();
    }

private:
    // Internal helper to write XML VTK StructuredGrid file (.vts)
    static void writeVTSFile(
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

        // XML header
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "  <StructuredGrid WholeExtent=\"0 " << (Nx-1) << " 0 " << (Ny-1) << " 0 0\">\n";
        file << "    <Piece Extent=\"0 " << (Nx-1) << " 0 " << (Ny-1) << " 0 0\">\n";
        
        // Points
        file << "      <Points>\n";
        file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                file << "          " << mesh.x(i) << " " << mesh.y(j) << " 0\n";
            }
        }
        file << "        </DataArray>\n";
        file << "      </Points>\n";
        
        // Point Data - specify default Scalars and Vectors for ParaView
        file << "      <PointData Scalars=\"velocity_magnitude\" Vectors=\"velocity\">\n";
        
        // Velocity magnitude
        file << "        <DataArray type=\"Float32\" Name=\"velocity_magnitude\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double mag = std::sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j));
                file << "          " << mag << "\n";
            }
        }
        file << "        </DataArray>\n";
        
        // Pressure
        file << "        <DataArray type=\"Float32\" Name=\"pressure\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                file << "          " << p(i,j) << "\n";
            }
        }
        file << "        </DataArray>\n";
        
        // Velocity vector
        file << "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                file << "          " << u(i,j) << " " << v(i,j) << " 0\n";
            }
        }
        file << "        </DataArray>\n";
        
        file << "      </PointData>\n";
        file << "    </Piece>\n";
        file << "  </StructuredGrid>\n";
        file << "</VTKFile>\n";

        if (!file.good()) {
            std::cerr << "[IO WARNING] Error during write to file: " << filename << std::endl;
        }
        file.close();
    }

    // Write VTS file with RGB color arrays based on ColorMap
    static void writeVTSFileWithColors(
        const std::string& filename,
        const core::Mesh& mesh,
        const core::Field& u,
        const core::Field& v,
        const core::Field& p,
        ColorScheme scheme = ColorScheme::BLUE_GREEN_RED)
    {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "[IO WARNING] Could not open file for writing: " << filename << std::endl;
            return;
        }

        int Nx = mesh.Nx();
        int Ny = mesh.Ny();

        // Compute velocity magnitude for coloring
        std::vector<double> vel_mag(Nx * Ny);
        double vmax = 0.0;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int idx = j * Nx + i;
                vel_mag[idx] = std::sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j));
                vmax = std::max(vmax, vel_mag[idx]);
            }
        }

        // Normalize and map to colors
        std::vector<RGB> colors(Nx * Ny);
        for (int idx = 0; idx < Nx * Ny; ++idx) {
            double normalized = (vmax > 1e-12) ? (vel_mag[idx] / vmax) : 0.0;
            colors[idx] = ColorMap::mapValue(normalized, scheme);
        }

        // XML header
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "  <StructuredGrid WholeExtent=\"0 " << (Nx-1) << " 0 " << (Ny-1) << " 0 0\">\n";
        file << "    <Piece Extent=\"0 " << (Nx-1) << " 0 " << (Ny-1) << " 0 0\">\n";
        
        // Points
        file << "      <Points>\n";
        file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                file << "          " << mesh.x(i) << " " << mesh.y(j) << " 0\n";
            }
        }
        file << "        </DataArray>\n";
        file << "      </Points>\n";
        
        // Point Data - keep analysis fields as defaults
        file << "      <PointData Scalars=\"velocity_magnitude\" Vectors=\"velocity\">\n";

        // Velocity magnitude (primary scalar for post-processing)
        file << "        <DataArray type=\"Float32\" Name=\"velocity_magnitude\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                file << "          " << vel_mag[j * Nx + i] << "\n";
            }
        }
        file << "        </DataArray>\n";
        
        // Pressure
        file << "        <DataArray type=\"Float32\" Name=\"pressure\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                file << "          " << p(i,j) << "\n";
            }
        }
        file << "        </DataArray>\n";
        
        // Velocity vector
        file << "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                file << "          " << u(i,j) << " " << v(i,j) << " 0\n";
            }
        }
        file << "        </DataArray>\n";

        // Optional helper RGB field for direct display (Blue->Green->Red)
        file << "        <DataArray type=\"UInt8\" Name=\"velocity_rgb\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int idx = j * Nx + i;
                file << "          "
                     << static_cast<int>(std::clamp(colors[idx].r, 0.0f, 1.0f) * 255.0f) << " "
                     << static_cast<int>(std::clamp(colors[idx].g, 0.0f, 1.0f) * 255.0f) << " "
                     << static_cast<int>(std::clamp(colors[idx].b, 0.0f, 1.0f) * 255.0f) << "\n";
            }
        }
        file << "        </DataArray>\n";
        
        file << "      </PointData>\n";
        file << "    </Piece>\n";
        file << "  </StructuredGrid>\n";
        file << "</VTKFile>\n";

        if (!file.good()) {
            std::cerr << "[IO WARNING] Error during write to file: " << filename << std::endl;
        }
        file.close();
    }

    // Write VTS file with a primary named scalar field (e.g. temperature) as default display
    static void writeVTSFileScalar(
        const std::string& filename,
        const core::Mesh& mesh,
        const core::Field& scalar,
        const std::string& scalarName,
        const core::Field& u,
        const core::Field& v)
    {
        std::ofstream file(filename);

        int Nx = mesh.Nx();
        int Ny = mesh.Ny();

        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "  <StructuredGrid WholeExtent=\"0 " << (Nx-1) << " 0 " << (Ny-1) << " 0 0\">\n";
        file << "    <Piece Extent=\"0 " << (Nx-1) << " 0 " << (Ny-1) << " 0 0\">\n";

        // Points
        file << "      <Points>\n";
        file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                file << "          " << mesh.x(i) << " " << mesh.y(j) << " 0\n";
        file << "        </DataArray>\n";
        file << "      </Points>\n";

        // PointData: set scalarName as default Scalars so ParaView shows it immediately
        file << "      <PointData Scalars=\"" << scalarName << "\" Vectors=\"velocity\">\n";

        // Primary scalar (e.g. temperature)
        file << "        <DataArray type=\"Float32\" Name=\"" << scalarName << "\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                file << "          " << scalar(i,j) << "\n";
        file << "        </DataArray>\n";

        // Velocity magnitude
        file << "        <DataArray type=\"Float32\" Name=\"velocity_magnitude\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                file << "          " << std::sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j)) << "\n";
        file << "        </DataArray>\n";

        // Velocity vector
        file << "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                file << "          " << u(i,j) << " " << v(i,j) << " 0\n";
        file << "        </DataArray>\n";

        file << "      </PointData>\n";
        file << "    </Piece>\n";
        file << "  </StructuredGrid>\n";
        file << "</VTKFile>\n";

        file.close();
    }

    // Detect existing VTS files in the output directory for resume
    void detectExistingFrames()
    {
        if (!std::filesystem::exists(outputDir)) {
            startFrame = 0;
            lastFrameFilename = "";
            return;
        }
        
        int maxFrame = -1;
        std::string pattern = baseName + "_";
        
        // Scan directory for VTS files matching our pattern (7-digit frame numbers)
        for (const auto& entry : std::filesystem::directory_iterator(outputDir)) {
            if (!entry.is_regular_file()) continue;
            
            std::string filename = entry.path().filename().string();
            
            // Check if filename matches pattern: {baseName}_NNNNNNN.vts (7 digits)
            if (filename.find(pattern) == 0 && filename.substr(filename.length() - 4) == ".vts") {
                // Extract frame number (7 digits before .vts)
                size_t underscorePos = filename.find_last_of('_');
                size_t dotPos = filename.find_last_of('.');
                if (underscorePos != std::string::npos && dotPos != std::string::npos) {
                    std::string frameStr = filename.substr(underscorePos + 1, dotPos - underscorePos - 1);
                    try {
                        int frameNum = std::stoi(frameStr);
                        if (frameNum > maxFrame) {
                            maxFrame = frameNum;
                            lastFrameFilename = outputDir + "/" + filename;
                        }
                    } catch (...) {
                        // Invalid frame number, skip
                    }
                }
            }
        }
        
        if (maxFrame >= 0) {
            startFrame = maxFrame + 1;  // Resume from next frame
            if (verbosity > 0) {
                std::cout << "  Detected existing frames 0-" << maxFrame << ", resuming from frame " << startFrame << std::endl;
            }
        } else {
            startFrame = 0;
            lastFrameFilename = "";
        }
    }
    
    // Load existing timestep entries from PVD file
    void loadExistingPVD()
    {
        std::string pvdFilename = (std::filesystem::path(pvdOutputDir) / (baseName + ".pvd")).string();
        if (!std::filesystem::exists(pvdFilename)) {
            return;
        }
        
        std::ifstream pvdFile(pvdFilename);
        if (!pvdFile.is_open()) {
            return;
        }
        
        std::string line;
        std::regex datasetRegex(R"REGEX(\s*<DataSet\s+timestep="([^"]+)"\s+file="([^"]+)"/>)REGEX");
        
        while (std::getline(pvdFile, line)) {
            std::smatch match;
            if (std::regex_search(line, match, datasetRegex)) {
                double time = std::stod(match[1].str());
                std::string filepath = match[2].str();
                
                // Extract filename from path (remove ./ prefix and directory)
                size_t lastSlash = filepath.find_last_of('/');
                std::string filename = (lastSlash != std::string::npos) 
                    ? filepath.substr(lastSlash + 1) 
                    : filepath;
                
                // Extract frame number from filename
                size_t underscorePos = filename.find_last_of('_');
                size_t dotPos = filename.find_last_of('.');
                if (underscorePos != std::string::npos && dotPos != std::string::npos) {
                    std::string frameStr = filename.substr(underscorePos + 1, dotPos - underscorePos - 1);
                    try {
                        int step = std::stoi(frameStr);
                        timesteps.push_back({step, time, filename});
                    } catch (...) {
                        // Skip invalid entries
                    }
                }
            }
        }
        
        pvdFile.close();
        
        if (verbosity > 0 && !timesteps.empty()) {
            std::cout << "  Loaded " << timesteps.size() << " existing timestep entries from PVD" << std::endl;
        }
    }
};

}
