/**
 * @file bench_gmp.cpp
 * @brief Performance comparison: bigint-cpp vs GMP (GNU Multiple Precision)
 *
 * This benchmark compares our BigInt implementation against GMP's mpz_t.
 * GMP is the gold standard for arbitrary precision arithmetic.
 *
 * Run with: ./bigint_gmp_benchmarks --benchmark_format=console
 *
 * To install GMP:
 *   - vcpkg: vcpkg install gmp:x64-windows
 *   - Ubuntu: sudo apt install libgmp-dev
 *   - macOS: brew install gmp
 */

#include <bigint/BigInt.hpp>
#include <gmp.h>

#include <benchmark/benchmark.h>

using namespace bigint;

// ============================================================================
// Test Numbers (same for both libraries)
// ============================================================================

static const char* NUM_SMALL = "12345678901234567890123456789012345678901234567890";
static const char* NUM_SMALL_B = "98765432109876543210987654321098765432109876543210";

static const char* NUM_LARGE =
    "123456789012345678901234567890123456789012345678901234567890"
    "123456789012345678901234567890123456789012345678901234567890"
    "123456789012345678901234567890123456789012345678901234567890";

static const char* NUM_LARGE_B =
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

static void BM_GMP_Add(benchmark::State& state) {
  mpz_t a, b, result;
  mpz_init_set_str(a, NUM_LARGE, 10);
  mpz_init_set_str(b, NUM_LARGE_B, 10);
  mpz_init(result);
  for (auto _ : state) {
    mpz_add(result, a, b);
    benchmark::DoNotOptimize(result);
  }
  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(result);
}
BENCHMARK(BM_GMP_Add);

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

static void BM_GMP_Multiply(benchmark::State& state) {
  mpz_t a, b, result;
  mpz_init_set_str(a, NUM_SMALL, 10);
  mpz_init_set_str(b, NUM_SMALL_B, 10);
  mpz_init(result);
  for (auto _ : state) {
    mpz_mul(result, a, b);
    benchmark::DoNotOptimize(result);
  }
  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(result);
}
BENCHMARK(BM_GMP_Multiply);

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

static void BM_GMP_Divide(benchmark::State& state) {
  mpz_t a, b, result;
  mpz_init_set_str(a, NUM_LARGE, 10);
  mpz_init_set_str(b, NUM_SMALL, 10);
  mpz_init(result);
  for (auto _ : state) {
    mpz_tdiv_q(result, a, b);
    benchmark::DoNotOptimize(result);
  }
  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(result);
}
BENCHMARK(BM_GMP_Divide);

// ============================================================================
// Modular Exponentiation Comparison (powmod)
// ============================================================================

static void BM_BigInt_PowMod(benchmark::State& state) {
  BigInt base(NUM_SMALL);
  BigInt exp("65537");
  BigInt mod(NUM_SMALL_B);
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod);

static void BM_GMP_PowMod(benchmark::State& state) {
  mpz_t base, exp, mod, result;
  mpz_init_set_str(base, NUM_SMALL, 10);
  mpz_init_set_str(exp, "65537", 10);
  mpz_init_set_str(mod, NUM_SMALL_B, 10);
  mpz_init(result);
  for (auto _ : state) {
    mpz_powm(result, base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
  mpz_clear(base);
  mpz_clear(exp);
  mpz_clear(mod);
  mpz_clear(result);
}
BENCHMARK(BM_GMP_PowMod);

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

static void BM_GMP_GCD(benchmark::State& state) {
  mpz_t a, b, result;
  mpz_init_set_str(a, NUM_LARGE, 10);
  mpz_init_set_str(b, NUM_LARGE_B, 10);
  mpz_init(result);
  for (auto _ : state) {
    mpz_gcd(result, a, b);
    benchmark::DoNotOptimize(result);
  }
  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(result);
}
BENCHMARK(BM_GMP_GCD);
