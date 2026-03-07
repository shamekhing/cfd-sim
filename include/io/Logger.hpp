#pragma once

#include <iostream>
#include <string>

namespace io {

class Logger {
public:
    inline static int verbosity = 1;

    static void setVerbosity(int v) { verbosity = v; }

    static void logInfo(const std::string& msg) {
        if (verbosity >= 1)
            std::cout << "[INFO] " << msg << std::endl;
    }

    static void logWarn(const std::string& msg) {
        if (verbosity >= 1)
            std::cerr << "[WARN] " << msg << std::endl;
    }

    static void logDebug(const std::string& msg) {
        if (verbosity >= 2)
            std::cout << "[DEBUG] " << msg << std::endl;
    }
};

}
