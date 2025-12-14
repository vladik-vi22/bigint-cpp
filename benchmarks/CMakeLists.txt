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
# Comparison Benchmark (vs Boost.Multiprecision - standalone mode)
# ============================================================================
option(BIGINT_BUILD_COMPARISON_BENCHMARKS "Build comparison benchmarks vs Boost" OFF)

if(BIGINT_BUILD_COMPARISON_BENCHMARKS)
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

    add_executable(bigint_comparison_benchmarks
        bench_comparison.cpp
    )

    target_compile_definitions(bigint_comparison_benchmarks
        PRIVATE
            BOOST_MP_STANDALONE
    )

    target_include_directories(bigint_comparison_benchmarks
        PRIVATE
            ${boost_config_SOURCE_DIR}/include
            ${boost_multiprecision_SOURCE_DIR}/include
    )

    target_link_libraries(bigint_comparison_benchmarks
        PRIVATE
            bigint::bigint
            benchmark::benchmark
            benchmark::benchmark_main
    )
endif()

