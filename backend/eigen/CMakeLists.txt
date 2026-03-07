cmake_minimum_required(VERSION 3.15)

# Eigen backend library
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

# Header-only interface library
add_library(cfd_backend_eigen INTERFACE)

target_include_directories(cfd_backend_eigen INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/backend/eigen>
)

# Link Eigen and OpenMP for parallel sparse operations
target_link_libraries(cfd_backend_eigen INTERFACE Eigen3::Eigen OpenMP::OpenMP_CXX)

# Make backend implementation available as cfd::backend::eigen
add_library(cfd::backend::eigen ALIAS cfd_backend_eigen)
