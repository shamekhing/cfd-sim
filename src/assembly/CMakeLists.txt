add_library(assembly STATIC
    GenericScalarAssembler.cpp
)

target_include_directories(assembly
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(assembly
    PUBLIC
    core
    linear
    numerics
)
