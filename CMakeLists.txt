cmake_minimum_required(VERSION 3.16)
project(CFD LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable testing for CTest
enable_testing()

option(CFD_SCALAR_FLOAT "Use float as global scalar type" OFF)
set(CFD_BACKEND "EIGEN" CACHE STRING "Linear algebra backend: EIGEN or CUDA")
set_property(CACHE CFD_BACKEND PROPERTY STRINGS EIGEN CUDA)

string(TOUPPER "${CFD_BACKEND}" CFD_BACKEND_UPPER)
if(NOT CFD_BACKEND_UPPER STREQUAL "EIGEN" AND NOT CFD_BACKEND_UPPER STREQUAL "CUDA")
    message(FATAL_ERROR "Invalid CFD_BACKEND='${CFD_BACKEND}'. Use EIGEN or CUDA.")
endif()

if(CFD_BACKEND_UPPER STREQUAL "CUDA")
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_ARCHITECTURES 80)  # Ampere (RTX 30 series, A100, etc.)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "CUDA backend enabled. CUDA Toolkit version: ${CUDAToolkit_VERSION}")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -march=native -ffast-math)
else()
    add_compile_options(-O0 -g)
endif()

add_compile_options(
    $<$<COMPILE_LANGUAGE:CXX>:-Wall>
    $<$<COMPILE_LANGUAGE:CXX>:-Wextra>
    $<$<COMPILE_LANGUAGE:CXX>:-Wpedantic>
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wall,-Wextra>
)

# Pass backend selection to compiler (0=Eigen, 1=CUDA)
if(CFD_BACKEND_UPPER STREQUAL "CUDA")
    add_compile_definitions(CFD_BACKEND=1)
    add_compile_definitions(CFD_ENABLE_CUDA=1)
else()
    add_compile_definitions(CFD_BACKEND=0)
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

# Build all modules
add_subdirectory(src/linear)
add_subdirectory(src/core)
add_subdirectory(src/numerics)
add_subdirectory(src/assembly)
add_subdirectory(src/physics)
add_subdirectory(src/solver)
add_subdirectory(src/config)
add_subdirectory(src/io)

# Build cases (examples)
add_subdirectory(cases/examples)

# Validation reference data
add_subdirectory(tests/validation)

file(GLOB_RECURSE SOURCES
    src/physics/*.cpp
)

if(CFD_BACKEND_UPPER STREQUAL "CUDA")
    file(GLOB_RECURSE CUDA_SOURCES
        linalg/backends/cuda/*.cu
    )
    set(SOURCES ${SOURCES} ${CUDA_SOURCES})
endif()

add_executable(cfd
    src/main.cpp
    ${SOURCES}
)

target_include_directories(cfd PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(cfd PRIVATE
    examples
    validation
    io
    config
    solver
    physics
    assembly
    numerics
    core
    linear
    Eigen3::Eigen
    Threads::Threads
)

target_compile_definitions(cfd PRIVATE
    $<$<BOOL:${CFD_SCALAR_FLOAT}>:CFD_SCALAR_FLOAT=1>
    $<$<NOT:$<BOOL:${CFD_SCALAR_FLOAT}>>:CFD_SCALAR_FLOAT=0>
    $<$<STREQUAL:${CFD_BACKEND_UPPER},CUDA>:CFD_BACKEND_CUDA=1>
    $<$<NOT:$<STREQUAL:${CFD_BACKEND_UPPER},CUDA>>:CFD_BACKEND_CUDA=0>
)

if(CFD_BACKEND_UPPER STREQUAL "CUDA")
    target_link_libraries(cfd PRIVATE CUDA::cusparse CUDA::cusolver CUDA::cudart)
    target_compile_options(cfd PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_80,code=sm_80>
    )
endif()

# Build tests
add_subdirectory(tests)
