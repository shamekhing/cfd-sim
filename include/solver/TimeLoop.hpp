#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>
#include <functional>
#include <utility>

namespace solver {

/**
 * General-purpose time stepping loop for transient simulations
 * 
 * Provides:
 *  - Time step iteration management
 *  - Progress tracking with elapsed/estimated time
 *  - Convergence monitoring (field residuals)
 *  - Statistics collection and reporting
 *  - Resume capability via start frame detection
 * 
 * Usage:
 * 
 *   TimeLoop loop(numTimeSteps, dt, verbosity);
 *   
 *   while (loop.isContinuing()) {
 *       int tstep = loop.getTimeStep();
 *       double time = loop.getTime();
 *       
 *       // Your solver call here
 *       solver(fields, time, dt);
 *       
 *       // Compute residual for convergence tracking
 *       double residual = computeResidual(fields);
 *       loop.recordResidual(residual);
 *       
 *       // Output frame (if needed)
 *       io.writeFrame(tstep, time, fields);
 *       
 *       loop.advance();
 *   }
 *   
 *   loop.printSummary();
 */
class TimeLoop {
private:
    int maxTimeSteps;
    double dt;
    double currentTime;
    int currentStep;
    int startStep;
    int verbosity;
    
    std::chrono::high_resolution_clock::time_point loopStart;
    std::chrono::high_resolution_clock::time_point loopCurrent;
    
    double maxResidual;
    double avgResidual;
    long residualCount;
    double maxFieldChange;

public:
    /**
     * Constructor
     * @param numSteps      Total number of time steps
     * @param timeStep      Time step size (dt)
     * @param verbosity     Output level: 0=silent, 1=progress bar, 2=detailed
     * @param initialStep   Starting step (for resume; default=0)
     */
    TimeLoop(int numSteps, double timeStep, int verbosity = 1, int initialStep = 0)
        : maxTimeSteps(numSteps),
          dt(timeStep),
          currentTime(initialStep * timeStep),
          currentStep(initialStep),
          startStep(initialStep),
          verbosity(verbosity),
          maxResidual(0.0),
          avgResidual(0.0),
          residualCount(0),
          maxFieldChange(0.0)
    {
        loopStart = std::chrono::high_resolution_clock::now();
        loopCurrent = loopStart;
    }

    // ============================================================
    // Main Loop Control
    // ============================================================

    /**
     * Check if the loop should continue
     */
    bool isContinuing() const {
        return currentStep < maxTimeSteps;
    }

    /**
     * Advance to next time step
     */
    void advance() {
        currentStep++;
        currentTime += dt;
        loopCurrent = std::chrono::high_resolution_clock::now();
        printProgress();
    }

    // ============================================================
    // Query Methods
    // ============================================================

    int getTimeStep() const { return currentStep; }
    double getTime() const { return currentTime; }
    double getTimeStepSize() const { return dt; }
    int getTotalSteps() const { return maxTimeSteps; }
    
    /**
     * Get elapsed wall-clock time in seconds
     */
    double getElapsedTime() const {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(now - loopStart).count();
    }

    /**
     * Get estimated time to completion in seconds
     */
    double getEstimatedRemaining() const {
        if (currentStep <= startStep) return 0.0;
        double elapsed = getElapsedTime();
        double avgTimePerStep = elapsed / (currentStep - startStep);
        int stepsRemaining = maxTimeSteps - currentStep;
        return avgTimePerStep * stepsRemaining;
    }

    // ============================================================
    // Convergence and Statistics
    // ============================================================

    /**
     * Record residual for convergence tracking
     * @param residual  Convergence metric (e.g., velocity divergence L2 norm)
     */
    void recordResidual(double residual) {
        maxResidual = std::max(maxResidual, residual);
        avgResidual += residual;
        residualCount++;
    }

    /**
     * Record maximum field change (for divergence-based convergence)
     */
    void recordFieldChange(double delta) {
        maxFieldChange = std::max(maxFieldChange, delta);
    }

    /**
     * Get statistics
     */
    struct Statistics {
        int stepsCompleted;
        double maxResidual;
        double avgResidual;
        double maxFieldChange;
        double wallTimeSeconds;
        double estRemainingSeconds;
    };

    Statistics getStatistics() const {
        Statistics stats;
        stats.stepsCompleted = currentStep - startStep;
        stats.maxResidual = maxResidual;
        stats.avgResidual = (residualCount > 0) ? avgResidual / residualCount : 0.0;
        stats.maxFieldChange = maxFieldChange;
        stats.wallTimeSeconds = getElapsedTime();
        stats.estRemainingSeconds = getEstimatedRemaining();
        return stats;
    }

    // ============================================================
    // Output and Diagnostics
    // ============================================================

    /**
     * Print iteration progress (called automatically by advance())
     */
    void printProgress() {
        if (verbosity < 1) return;
        if (currentStep == startStep && verbosity < 2) return;

        auto now = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(now - loopStart).count();
        double avgTime = (currentStep > startStep) 
            ? elapsed / (currentStep - startStep) 
            : 0.0;
        double remaining = avgTime * (maxTimeSteps - currentStep);

        int progress = (currentStep * 100) / (maxTimeSteps - 1);
        std::cout << "\r  Timestep " << currentStep << "/" << maxTimeSteps 
                  << " [" << std::string(progress / 5, '#') 
                  << std::string(20 - progress / 5, '-') << "] " 
                  << progress << "% | "
                  << std::fixed << std::setprecision(1)
                  << "T:" << elapsed << "s | ETA:" << remaining << "s";
        
        if (residualCount > 0 && verbosity >= 2) {
            std::cout << " | Resid:" << std::scientific << std::setprecision(2) 
                      << (maxResidual > 0 ? maxResidual : 0.0) << std::defaultfloat;
        }
        
        std::cout << std::flush;
        
        if (currentStep == maxTimeSteps - 1) {
            std::cout << std::endl;
        }
    }

    /**
     * Print summary statistics at end of simulation
     */
    void printSummary() {
        if (verbosity < 1) return;

        double elapsed = getElapsedTime();
        double avgTimePerStep = (currentStep > startStep) 
            ? elapsed / (currentStep - startStep) 
            : 0.0;

        std::cout << "\n  ===== TimeLoop Summary =====" << std::endl;
        std::cout << "  Steps: " << (currentStep - startStep) << " / " << maxTimeSteps 
                  << " completed" << std::endl;
        std::cout << "  Total time: " << std::fixed << std::setprecision(2) 
                  << elapsed << "s" << std::endl;
        std::cout << "  Avg/step: " << avgTimePerStep << "s" << std::endl;
        
        if (residualCount > 0) {
            std::cout << "  Max residual: " << std::scientific << std::setprecision(2) 
                      << maxResidual << std::defaultfloat << std::endl;
            std::cout << "  Avg residual: " << std::scientific << std::setprecision(2) 
                      << (avgResidual / residualCount) << std::defaultfloat << std::endl;
        }
        
        if (maxFieldChange > 0) {
            std::cout << "  Max field change: " << std::scientific << std::setprecision(2) 
                      << maxFieldChange << std::defaultfloat << std::endl;
        }
        
        std::cout << "  ============================\n" << std::endl;
    }

};

}  // namespace solver