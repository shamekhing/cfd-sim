cmake_minimum_required(VERSION 3.15)

# Backend layer - contains ALL vendor/architecture-specific code
# Application code must use backend/BackendDispatcher.hpp to route to active backend

if(NOT CFD_BACKEND)
    set(CFD_BACKEND "EIGEN" CACHE STRING "Active backend: EIGEN or CUDA")
endif()

if(CFD_BACKEND STREQUAL "CUDA")
    message(STATUS "Backend: CUDA")
    add_subdirectory(cuda)
else()
    message(STATUS "Backend: Eigen (default)")
    add_subdirectory(eigen)
endif()
