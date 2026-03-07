add_library(physics STATIC
    ScalarTransport.cpp
)

target_include_directories(physics
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(physics
    PUBLIC
    assembly
    numerics
    core
    linear
)
