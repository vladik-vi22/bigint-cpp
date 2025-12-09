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

