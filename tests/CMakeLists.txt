cmake_minimum_required(VERSION 3.15)

enable_testing()

# Unit tests (MMS verification tests)
add_subdirectory(unit)

# Note: Validation tests are defined in tests/validation/CMakeLists.txt
# Examples (simulation cases) are in cases/examples/
