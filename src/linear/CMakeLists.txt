cmake_minimum_required(VERSION 3.15)

set(LINEAR_SOURCES
    Vector.cpp
    Matrix.cpp
    Solver.cpp
    BiCGSTAB.cpp
    CG.cpp
    LU.cpp
    AdaptiveSolver.cpp
)

# Add CUDA backend if enabled
if(CFD_BACKEND_UPPER STREQUAL "CUDA")
    list(APPEND LINEAR_SOURCES ../../backend/cuda/src/CudaBackend.cu)
    enable_language(CUDA)
endif()

add_library(linear STATIC ${LINEAR_SOURCES})

# Include directories for angle-bracket includes like <linear/Vector.hpp>
target_include_directories(linear PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../backend/cuda/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../backend
)

# Eigen dependency
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
target_link_libraries(linear PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX)

# CUDA backend (optional)
if(CFD_BACKEND_UPPER STREQUAL "CUDA")
    target_link_libraries(linear PUBLIC 
        CUDA::cublas
        CUDA::cusparse
        CUDA::cudart
    )
    target_compile_definitions(linear PUBLIC CFD_BACKEND_CUDA=1)
else()
    target_compile_definitions(linear PUBLIC CFD_BACKEND_CUDA=0)
endif()
