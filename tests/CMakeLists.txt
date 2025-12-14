# ============================================================================
# GoogleTest Dependency
# ============================================================================
include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# ============================================================================
# Test Executable
# ============================================================================
add_executable(bigint_tests
    test_bigint_basic.cpp
    test_bigint_arithmetic.cpp
    test_bigint_bitwise.cpp
    test_bigint_properties.cpp
)

target_link_libraries(bigint_tests
    PRIVATE
        bigint::bigint
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(bigint_tests)

