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

// ~8640 bits (270 words) - triggers Toom-Cook 3-way multiplication
static std::string generateToom3String() {
  std::string result;
  for (int i = 0; i < 32; ++i) {
    result += "12345678901234567890123456789012345678901234567890123456789012345678901234567890";
  }
  return result;
}

static std::string generateToom3StringB() {
  std::string result;
  for (int i = 0; i < 32; ++i) {
    result += "98765432109876543210987654321098765432109876543210987654321098765432109876543210";
  }
  return result;
}

static const std::string NUM_TOOM3 = generateToom3String();
static const std::string NUM_TOOM3_B = generateToom3StringB();

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

static void BM_BigInt_Multiply_Toom3(benchmark::State& state) {
  BigInt a(NUM_TOOM3, 10);
  BigInt b(NUM_TOOM3_B, 10);
  for (auto _ : state) {
    BigInt result = a * b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Multiply_Toom3);

static void BM_Boost_Multiply_Toom3(benchmark::State& state) {
  cpp_int a(NUM_TOOM3);
  cpp_int b(NUM_TOOM3_B);
  for (auto _ : state) {
    cpp_int result = a * b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_Boost_Multiply_Toom3);

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
// Modular Exponentiation Comparison (powmod) - RSA Public Key (small exponent)
// ============================================================================

static void BM_BigInt_PowMod(benchmark::State& state) {
  BigInt base(NUM_SMALL);
  BigInt exp("65537");  // RSA public exponent (17 bits)
  BigInt mod(NUM_SMALL_B);  // Even modulus - uses standard algorithm
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
// Modular Exponentiation - RSA Private Key Simulation (large exponent)
// ============================================================================

// Large odd modulus for Montgomery algorithm testing (~1024 bits)
static const std::string NUM_LARGE_ODD =
    "1234567890123456789012345678901234567890123456789012345678901234567890"
    "1234567890123456789012345678901234567890123456789012345678901234567890"
    "1234567890123456789012345678901234567890123456789012345678901234567890"
    "1234567890123456789012345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567891";

static void BM_BigInt_PowMod_LargeExp(benchmark::State& state) {
  BigInt base(NUM_SMALL);
  BigInt exp(NUM_SMALL);  // Large exponent (~200 bits) - simulates private key
  BigInt mod(NUM_LARGE_ODD);  // Large odd modulus - uses Montgomery CIOS
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_LargeExp);

static void BM_Boost_PowMod_LargeExp(benchmark::State& state) {
  cpp_int base(NUM_SMALL);
  cpp_int exp(NUM_SMALL);  // Large exponent (~200 bits)
  cpp_int mod(NUM_LARGE_ODD);
  for (auto _ : state) {
    cpp_int result = powm(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_Boost_PowMod_LargeExp);

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

// ============================================================================
// Divmod Comparison (quotient + remainder in one operation)
// ============================================================================

static void BM_BigInt_Divmod(benchmark::State& state) {
  BigInt a(NUM_LARGE);
  BigInt b(NUM_SMALL);
  for (auto _ : state) {
    auto [q, r] = divmod(a, b);
    benchmark::DoNotOptimize(q);
    benchmark::DoNotOptimize(r);
  }
}
BENCHMARK(BM_BigInt_Divmod);

static void BM_Boost_Divmod(benchmark::State& state) {
  cpp_int a(NUM_LARGE);
  cpp_int b(NUM_SMALL);
  for (auto _ : state) {
    cpp_int q, r;
    boost::multiprecision::divide_qr(a, b, q, r);
    benchmark::DoNotOptimize(q);
    benchmark::DoNotOptimize(r);
  }
}
BENCHMARK(BM_Boost_Divmod);

// ============================================================================
// PopCount Comparison (count set bits)
// ============================================================================

static void BM_BigInt_PopCount(benchmark::State& state) {
  BigInt a(NUM_LARGE);
  for (auto _ : state) {
    size_t result = a.popCount();
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PopCount);

static void BM_Boost_PopCount(benchmark::State& state) {
  cpp_int a(NUM_LARGE);
  for (auto _ : state) {
    // Boost doesn't have popcount - count bits manually via bit_test
    unsigned count = 0;
    unsigned bits = msb(a) + 1;
    for (unsigned i = 0; i < bits; ++i) {
      if (bit_test(a, i)) ++count;
    }
    benchmark::DoNotOptimize(count);
  }
}
BENCHMARK(BM_Boost_PopCount);

// ============================================================================
// TestBit Comparison
// ============================================================================

static void BM_BigInt_TestBit(benchmark::State& state) {
  BigInt a(NUM_LARGE);
  for (auto _ : state) {
    bool result = a.testBit(100);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_TestBit);

static void BM_Boost_TestBit(benchmark::State& state) {
  cpp_int a(NUM_LARGE);
  for (auto _ : state) {
    bool result = bit_test(a, 100);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_Boost_TestBit);

// ============================================================================
// SetBit Comparison
// ============================================================================

static void BM_BigInt_SetBit(benchmark::State& state) {
  BigInt a(NUM_SMALL);
  for (auto _ : state) {
    a.setBit(100);
    a.clearBit(100);  // Reset for next iteration
    benchmark::DoNotOptimize(a);
  }
}
BENCHMARK(BM_BigInt_SetBit);

static void BM_Boost_SetBit(benchmark::State& state) {
  cpp_int a(NUM_SMALL);
  for (auto _ : state) {
    bit_set(a, 100);
    bit_unset(a, 100);  // Reset for next iteration
    benchmark::DoNotOptimize(a);
  }
}
BENCHMARK(BM_Boost_SetBit);

// ============================================================================
// FlipBit Comparison
// ============================================================================

static void BM_BigInt_FlipBit(benchmark::State& state) {
  BigInt a(NUM_SMALL);
  for (auto _ : state) {
    a.flipBit(50);
    benchmark::DoNotOptimize(a);
  }
}
BENCHMARK(BM_BigInt_FlipBit);

static void BM_Boost_FlipBit(benchmark::State& state) {
  cpp_int a(NUM_SMALL);
  for (auto _ : state) {
    bit_flip(a, 50);
    benchmark::DoNotOptimize(a);
  }
}
BENCHMARK(BM_Boost_FlipBit);

BENCHMARK_MAIN();

