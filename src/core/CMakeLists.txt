add_library(core
    Mesh.cpp
    Field.cpp
    BoundaryCondition.cpp
)

target_include_directories(core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../..
)

target_link_libraries(core PUBLIC linear)

target_compile_features(core PUBLIC cxx_std_17)
