/**
 * @file bench_comparison.cpp
 * @brief Performance comparison: bigint-cpp vs Boost.Multiprecision
 *
 * This benchmark compares our BigInt implementation against Boost's cpp_int.
 * Run with: ./bigint_comparison_benchmarks --benchmark_format=console
 *
 * Expected results: Boost is typically faster (highly optimized, assembly),
 * but bigint-cpp should be competitive for educational/prototyping use.
 */

#include <bigint/BigInt.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <benchmark/benchmark.h>

using namespace bigint;
using boost::multiprecision::cpp_int;

// ============================================================================
// Test Numbers (same for both libraries)
// ============================================================================

// ~200 bits
static const std::string NUM_SMALL = "12345678901234567890123456789012345678901234567890";
static const std::string NUM_SMALL_B = "98765432109876543210987654321098765432109876543210";

// ~600 bits
static const std::string NUM_LARGE =
    "123456789012345678901234567890123456789012345678901234567890"
    "123456789012345678901234567890123456789012345678901234567890"
    "123456789012345678901234567890123456789012345678901234567890";

static const std::string NUM_LARGE_B =
    "987654321098765432109876543210987654321098765432109876543210"
    "987654321098765432109876543210987654321098765432109876543210"
    "987654321098765432109876543210987654321098765432109876543210";

// ============================================================================
// Addition Comparison
// ============================================================================

static void BM_BigInt_Add(benchmark::State& state) {
  BigInt a(NUM_LARGE);
  BigInt b(NUM_LARGE_B);
  for (auto _ : state) {
    BigInt result = a + b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Add);

static void BM_Boost_Add(benchmark::State& state) {
  cpp_int a(NUM_LARGE);
  cpp_int b(NUM_LARGE_B);
  for (auto _ : state) {
    cpp_int result = a + b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_Boost_Add);

// ============================================================================
// Multiplication Comparison
// ============================================================================

static void BM_BigInt_Multiply(benchmark::State& state) {
  BigInt a(NUM_SMALL);
  BigInt b(NUM_SMALL_B);
  for (auto _ : state) {
    BigInt result = a * b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Multiply);

static void BM_Boost_Multiply(benchmark::State& state) {
  cpp_int a(NUM_SMALL);
  cpp_int b(NUM_SMALL_B);
  for (auto _ : state) {
    cpp_int result = a * b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_Boost_Multiply);

// ============================================================================
// Division Comparison
// ============================================================================

static void BM_BigInt_Divide(benchmark::State& state) {
  BigInt a(NUM_LARGE);
  BigInt b(NUM_SMALL);
  for (auto _ : state) {
    BigInt result = a / b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Divide);

static void BM_Boost_Divide(benchmark::State& state) {
  cpp_int a(NUM_LARGE);
  cpp_int b(NUM_SMALL);
  for (auto _ : state) {
    cpp_int result = a / b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_Boost_Divide);

// ============================================================================
// Modular Exponentiation Comparison (powmod)
// ============================================================================

static void BM_BigInt_PowMod(benchmark::State& state) {
  BigInt base(NUM_SMALL);
  BigInt exp("65537");  // RSA public exponent
  BigInt mod(NUM_SMALL_B);
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod);

static void BM_Boost_PowMod(benchmark::State& state) {
  cpp_int base(NUM_SMALL);
  cpp_int exp("65537");
  cpp_int mod(NUM_SMALL_B);
  for (auto _ : state) {
    cpp_int result = powm(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_Boost_PowMod);

// ============================================================================
// GCD Comparison
// ============================================================================

static void BM_BigInt_GCD(benchmark::State& state) {
  BigInt a(NUM_LARGE);
  BigInt b(NUM_LARGE_B);
  for (auto _ : state) {
    BigInt result = gcd(a, b);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_GCD);

static void BM_Boost_GCD(benchmark::State& state) {
  cpp_int a(NUM_LARGE);
  cpp_int b(NUM_LARGE_B);
  for (auto _ : state) {
    cpp_int result = boost::multiprecision::gcd(a, b);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_Boost_GCD);

BENCHMARK_MAIN();

