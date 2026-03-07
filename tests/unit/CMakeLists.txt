# Unit tests with known analytical solutions (MMS verification)

# Test 1: Poisson solver (Laplacian operator)
add_executable(test_poisson test_poisson.cpp)
target_include_directories(test_poisson PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)
target_link_libraries(test_poisson PRIVATE
    core
    assembly
    linear
)

# Test 2: Scalar transport with MMS
add_executable(test_scalar_mms test_scalar_mms.cpp)
target_include_directories(test_scalar_mms PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)
target_link_libraries(test_scalar_mms PRIVATE
    core
    assembly
    linear
)

# Register with CTest
add_test(NAME PoissonMMS COMMAND test_poisson)
add_test(NAME ScalarTransportMMS COMMAND test_scalar_mms)

set_tests_properties(PoissonMMS PROPERTIES
    TIMEOUT 60
    LABELS "unit;mms;poisson"
)

set_tests_properties(ScalarTransportMMS PROPERTIES
    TIMEOUT 120
    LABELS "unit;mms;scalar"
)
