add_library(examples INTERFACE)

target_include_directories(examples INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)

# Examples are header-only files that can be included by main.cpp
# They provide runLidDrivenCavity(), runHighReynoldsCavity(), etc.
