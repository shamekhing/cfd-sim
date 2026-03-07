cmake_minimum_required(VERSION 3.15)

# CUDA backend library (stub for now)

# Header-only interface library
add_library(cfd_backend_cuda INTERFACE)

target_include_directories(cfd_backend_cuda INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/backend/cuda>
)

# Optional CUDA support
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    target_link_libraries(cfd_backend_cuda INTERFACE CUDA::cudart CUDA::cublas CUDA::cusparse)
    
    # Add CUDA implementation sources when available
    target_sources(cfd_backend_cuda INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CudaBackend.cu>
    )
endif()

# Make backend available as cfd::backend::cuda
add_library(cfd::backend::cuda ALIAS cfd_backend_cuda)
