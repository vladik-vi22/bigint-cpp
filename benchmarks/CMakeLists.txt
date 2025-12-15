# ============================================================================
# Google Benchmark Dependency
# ============================================================================
include(FetchContent)

FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.3
)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(benchmark)

# ============================================================================
# Benchmark Executable
# ============================================================================
add_executable(bigint_benchmarks
    bench_bigint.cpp
)

target_link_libraries(bigint_benchmarks
    PRIVATE
        bigint::bigint
        benchmark::benchmark
        benchmark::benchmark_main
)

# ============================================================================
# Boost Comparison Benchmark (vs Boost.Multiprecision - standalone mode)
# ============================================================================
option(BIGINT_BUILD_BOOST_BENCHMARKS "Build comparison benchmarks vs Boost.Multiprecision" OFF)

if(BIGINT_BUILD_BOOST_BENCHMARKS)
    # Fetch Boost.Multiprecision and its minimal dependencies for standalone mode
    FetchContent_Declare(
        boost_config
        GIT_REPOSITORY https://github.com/boostorg/config.git
        GIT_TAG boost-1.84.0
    )
    FetchContent_Declare(
        boost_multiprecision
        GIT_REPOSITORY https://github.com/boostorg/multiprecision.git
        GIT_TAG boost-1.84.0
    )
    FetchContent_MakeAvailable(boost_config boost_multiprecision)

    add_executable(bigint_boost_benchmarks
        bench_boost.cpp
    )

    target_compile_definitions(bigint_boost_benchmarks
        PRIVATE
            BOOST_MP_STANDALONE
    )

    target_include_directories(bigint_boost_benchmarks
        PRIVATE
            ${boost_config_SOURCE_DIR}/include
            ${boost_multiprecision_SOURCE_DIR}/include
    )

    target_link_libraries(bigint_boost_benchmarks
        PRIVATE
            bigint::bigint
            benchmark::benchmark
            benchmark::benchmark_main
    )
endif()


# ============================================================================
# GMP Comparison Benchmark (requires GMP to be installed)
# ============================================================================
option(BIGINT_BUILD_GMP_BENCHMARKS "Build comparison benchmarks vs GMP (requires GMP)" OFF)

if(BIGINT_BUILD_GMP_BENCHMARKS)
    # Try to find GMP - user must have it installed (e.g., via vcpkg, MSYS2, or system)
    find_library(GMP_LIBRARY NAMES gmp libgmp mpir HINTS ENV GMP_ROOT)
    find_path(GMP_INCLUDE_DIR NAMES gmp.h HINTS ENV GMP_ROOT)

    if(GMP_LIBRARY AND GMP_INCLUDE_DIR)
        message(STATUS "Found GMP: ${GMP_LIBRARY}")
        message(STATUS "GMP include: ${GMP_INCLUDE_DIR}")

        add_executable(bigint_gmp_benchmarks
            bench_gmp.cpp
        )

        target_include_directories(bigint_gmp_benchmarks
            PRIVATE
                ${GMP_INCLUDE_DIR}
        )

        target_link_libraries(bigint_gmp_benchmarks
            PRIVATE
                bigint::bigint
                ${GMP_LIBRARY}
                benchmark::benchmark
                benchmark::benchmark_main
        )
    else()
        message(WARNING "GMP not found. Skipping GMP benchmarks.")
        message(STATUS "To install GMP:")
        message(STATUS "  - vcpkg: vcpkg install gmp:x64-windows")
        message(STATUS "  - Ubuntu: sudo apt install libgmp-dev")
        message(STATUS "  - macOS: brew install gmp")
    endif()
endif()
