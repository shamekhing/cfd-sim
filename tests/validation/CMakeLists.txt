# Validation reference data library
add_library(validation STATIC
    GhiaReference.cpp
)

target_include_directories(validation
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(validation
    PUBLIC
    core
    linear
)

# Validation test executables
# Test 1: Ghia Re=100 benchmark
add_executable(test_cavity_re100 test_cavity_re100.cpp)
target_include_directories(test_cavity_re100 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(test_cavity_re100 PRIVATE
    validation
    core
    solver
    numerics
    linear
    io
)

# Test 2: Ghia Re=1000 benchmark
add_executable(test_cavity_re1000 test_cavity_re1000.cpp)
target_include_directories(test_cavity_re1000 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(test_cavity_re1000 PRIVATE
    validation
    core
    solver
    numerics
    linear
    io
)

# Register with CTest
add_test(NAME CavityRe100_Ghia1982 COMMAND test_cavity_re100)
add_test(NAME CavityRe1000_Ghia1982 COMMAND test_cavity_re1000)

# Set test properties
set_tests_properties(CavityRe100_Ghia1982 PROPERTIES
    TIMEOUT 300  # 5 minutes
    LABELS "validation;benchmark;cavity"
)

set_tests_properties(CavityRe1000_Ghia1982 PROPERTIES
    TIMEOUT 600  # 10 minutes (Re=1000 takes longer)
    LABELS "validation;benchmark;cavity"
)
