add_library(solver INTERFACE)

target_include_directories(solver
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)

target_link_libraries(solver
    INTERFACE
    physics
    assembly
    numerics
    core
    linear
)
