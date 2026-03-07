#pragma once

// Backend dispatcher - routes to active backend at compile time
// Application code should ONLY include this file, NOT individual backend headers
// The dispatcher handles all backend selection logic

#ifndef CFD_BACKEND
#define CFD_BACKEND 0  // 0=Eigen, 1=CUDA
#endif

// Try to use requested backend
#if CFD_BACKEND == 1  // CUDA requested
    #ifdef CFD_ENABLE_CUDA
        #include "cuda/include/CudaBackend.hpp"
        namespace BackendImpl = cfd::backend::cuda;
    #else
        #error "CUDA backend requested but CFD_ENABLE_CUDA not defined. Set CFD_ENABLE_CUDA or use CFD_BACKEND=0"
    #endif
#else  // Default to Eigen
    #include "eigen/include/EigenBackend.hpp"
    namespace BackendImpl = cfd::backend::eigen;
#endif


