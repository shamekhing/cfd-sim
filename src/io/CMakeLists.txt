add_library(io INTERFACE)

target_include_directories(io
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)

target_link_libraries(io
    INTERFACE
    core
    linear
)
