#include <string>
#include <iostream>
#include <cstdlib>
#include "include/linear/Vector.hpp"
#include "include/config/SolverConfig.hpp"
#include "include/config/CaseConfigs.hpp"
#include "cases/examples/LidDrivenCavity.hpp"
#include "cases/examples/HighReynoldsCavity.hpp"
#include "cases/examples/ChannelFlow.hpp"
#include "cases/examples/EnergyTransport.hpp"
#include "cases/examples/CylinderFlow.hpp"
#include "cases/main/CylinderMain.hpp"

int main(int argc, char* argv[])
{
    // Print backend information
#if CFD_BACKEND == 1
    std::cout << "=== CFD Engine ===\nBackend: CUDA (GPU acceleration)\n" << std::endl;
#else
    std::cout << "=== CFD Engine ===\nBackend: Eigen (CPU with OpenMP)\n" << std::endl;
#endif

    linear::Scalar Re = 100.0;
    bool re_set = false;
    int Nx = 64;
    int Ny = 64;
    bool nx_set = false;
    bool ny_set = false;
    std::string caseType = "lid_cavity";  // default case
    std::string animationName = "output";  // default animation/output name
    std::string outputDir = "build/output";  // default output directory
    int numTimeSteps = 200;  // default timesteps for animation
    bool steps_set = false;
    bool enable_turbulence = false;  // turbulence modeling flag
    bool validate_mode = false;
    bool animate = true;  // global animation flag - enabled by default
    
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--re" && i + 1 < argc)
        {
            Re = std::atof(argv[++i]);
            re_set = true;
        }
        else if (arg == "--nx" && i + 1 < argc)
        {
            Nx = std::atoi(argv[++i]);
            nx_set = true;
        }
        else if (arg == "--ny" && i + 1 < argc)
        {
            Ny = std::atoi(argv[++i]);
            ny_set = true;
        }
        else if (arg == "--case" && i + 1 < argc)
        {
            caseType = argv[++i];
        }
        else if (arg == "--animate" || arg == "--anim")
        {
            animate = true;
        }
        else if ((arg == "--name" || arg == "-n") && i + 1 < argc)
        {
            animationName = argv[++i];
        }
        else if ((arg == "--output" || arg == "-o") && i + 1 < argc)
        {
            outputDir = argv[++i];
        }
        else if ((arg == "--numTimeSteps" || arg == "--steps") && i + 1 < argc)
        {
            numTimeSteps = std::atoi(argv[++i]);
            steps_set = true;
        }
        else if (arg == "--turbulence" || arg == "--turb")
        {
            enable_turbulence = true;
        }
        else if (arg == "--validate")
        {
            validate_mode = true;
        }
        else if (arg == "--help")
        {
            std::cout << "CFD Engine\n"
                      << "Usage: cfd [options]\n"
                      << "Options:\n"
                      << "  --case <name>       Case to run (default: lid_cavity)\n"
                      << "                      - lid_cavity: Benchmark cavity flow\n"
                      << "                      - high_re_cavity: Turbulent cavity (Re>500)\n"
                      << "                      - channel_flow: Periodic channel flow\n"
                      << "                      - energy_transport: Energy MMS validation\n"
                      << "                      - cylinder_flow: Vortex shedding (circle in flow)\n"
                      << "                      - cylinder_main: Flagship cylinder (dead-centre, Williamson)\n"
                      << "  --re <value>        Reynolds number (default: case-specific)\n"
                      << "  --nx <value>        Grid cells in x direction\n"
                      << "  --ny <value>        Grid cells in y direction\n"
                      << "  --animate           Enable time-stepping animation output\n"
                      << "  --name <name>       Output folder/file base name (default: output)\n"
                      << "  --output <path>     Output directory (default: build/output)\n"
                      << "  --numTimeSteps <n>  Number of time steps for animation (default: 200)\n"
                      << "  --validate          Use strict validation parameters\n"
                      << "  --turbulence        Enable k-epsilon turbulence model\n"
                      << "  --help              Show this help message\n";
            return 0;
        }
    }

    if (caseType == "energy_transport") {
        config::EnergyTransportConfig cfg;
        cfg.Nx = nx_set ? Nx : 64;
        cfg.Ny = ny_set ? Ny : 64;
        if (steps_set) cfg.numTimeSteps = numTimeSteps;
        cfg.useAnimation = animate;
        cases::runEnergyTransport(cfg, outputDir, animationName);
    } else if (caseType == "cylinder_flow") {
        config::CylinderFlowConfig cfg;
        cfg.Re = Re;
        cfg.Nx = nx_set ? Nx : cfg.Nx;
        cfg.Ny = ny_set ? Ny : cfg.Ny;
        if (steps_set) cfg.numTimeSteps = numTimeSteps;
        cfg.useAnimation = animate;
        cases::runCylinderFlow(cfg, outputDir, animationName);
    } else if (caseType == "cylinder_main") {
        config::CylinderFlowConfig cfg;
        if (re_set) cfg.Re = Re;
        cfg.Nx = nx_set ? Nx : cfg.Nx;
        cfg.Ny = ny_set ? Ny : cfg.Ny;
        if (steps_set) cfg.numTimeSteps = numTimeSteps;
        cfg.useAnimation = animate;
        cases::runCylinderMain(cfg, outputDir, animationName.empty() ? "cylinder" : animationName);
    } else if (caseType == "high_re_cavity") {
        config::HighReynoldsCavityConfig cfg;
        if (re_set) cfg.Re = Re;
        cfg.Nx = nx_set ? Nx : 64;
        cfg.Ny = ny_set ? Ny : 64;
        if (steps_set) cfg.numTimeSteps = numTimeSteps;
        cfg.useAnimation = animate;
        cases::runHighReynoldsCavity(cfg, outputDir, animationName);
    } else if (caseType == "channel_flow") {
        config::ChannelFlowConfig cfg;
        if (re_set) cfg.Re_tau = Re;
        cfg.Nx = nx_set ? Nx : 128;
        cfg.Ny = ny_set ? Ny : 64;
        if (steps_set) cfg.numTimeSteps = numTimeSteps;
        cfg.useAnimation = animate;
        cases::runChannelFlow(cfg, outputDir, animationName);
    } else if (caseType == "lid_cavity") {
        config::LidDrivenCavityConfig cfg;
        cfg.Re = Re;
        cfg.Nx = nx_set ? Nx : 64;
        cfg.Ny = ny_set ? Ny : 64;
        if (steps_set) cfg.numTimeSteps = numTimeSteps;
        cfg.useAnimation = animate;
        cfg.enable_turbulence = enable_turbulence;
        cfg.validationMode = validate_mode;
        cases::runLidDrivenCavity(cfg, outputDir, animationName);
    } else {
        std::cerr << "Unknown case: " << caseType << "\n";
        return 1;
    }

    return 0;
}
