add_library(numerics
    ConvectionSchemes.cpp
)

target_include_directories(numerics PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../..
)

target_link_libraries(numerics PUBLIC core linear)

target_compile_features(numerics PUBLIC cxx_std_17)
