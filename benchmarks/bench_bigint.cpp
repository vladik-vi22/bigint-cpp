#include <bigint/BigInt.hpp>

#include <benchmark/benchmark.h>

using namespace bigint;

// ============================================================================
// Addition Benchmarks
// ============================================================================

static void BM_BigInt_Addition_Small(benchmark::State& state) {
  BigInt a("12345678901234567890");
  BigInt b("98765432109876543210");
  for (auto _ : state) {
    BigInt result = a + b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Addition_Small);

static void BM_BigInt_Addition_Large(benchmark::State& state) {
  BigInt a("123456789012345678901234567890123456789012345678901234567890");
  BigInt b("987654321098765432109876543210987654321098765432109876543210");
  for (auto _ : state) {
    BigInt result = a + b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Addition_Large);

// ============================================================================
// Multiplication Benchmarks
// ============================================================================

static void BM_BigInt_Multiplication_Small(benchmark::State& state) {
  BigInt a("12345678901234567890");
  BigInt b("98765432109876543210");
  for (auto _ : state) {
    BigInt result = a * b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Multiplication_Small);

static void BM_BigInt_Multiplication_Large(benchmark::State& state) {
  BigInt a("123456789012345678901234567890123456789012345678901234567890");
  BigInt b("987654321098765432109876543210987654321098765432109876543210");
  for (auto _ : state) {
    BigInt result = a * b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Multiplication_Large);

// ============================================================================
// Division Benchmarks
// ============================================================================

static void BM_BigInt_Division_Small(benchmark::State& state) {
  BigInt a("98765432109876543210");
  BigInt b("12345678901234567890");
  for (auto _ : state) {
    BigInt result = a / b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Division_Small);

static void BM_BigInt_Division_Large(benchmark::State& state) {
  BigInt a("987654321098765432109876543210987654321098765432109876543210");
  BigInt b("123456789012345678901234567890");
  for (auto _ : state) {
    BigInt result = a / b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Division_Large);

// ============================================================================
// Modular Exponentiation Benchmarks (powmod)
// ============================================================================

static void BM_BigInt_PowMod_Small(benchmark::State& state) {
  BigInt base("12345");
  BigInt exp("100");
  BigInt mod("1000000007");
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_Small);

static void BM_BigInt_PowMod_Medium(benchmark::State& state) {
  BigInt base("123456789");
  BigInt exp("1000");
  BigInt mod("1000000007");
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_Medium);

static void BM_BigInt_PowMod_Large(benchmark::State& state) {
  BigInt base("123456789012345678901234567890");
  BigInt exp("10000");
  BigInt mod("1000000000000000000000000000057");
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_Large);

// ============================================================================
// GCD Benchmarks
// ============================================================================

static void BM_BigInt_GCD_Small(benchmark::State& state) {
  BigInt a("12345678901234567890");
  BigInt b("98765432109876543210");
  for (auto _ : state) {
    BigInt result = gcd(a, b);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_GCD_Small);

static void BM_BigInt_GCD_Large(benchmark::State& state) {
  BigInt a("123456789012345678901234567890123456789012345678901234567890");
  BigInt b("987654321098765432109876543210987654321098765432109876543210");
  for (auto _ : state) {
    BigInt result = gcd(a, b);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_GCD_Large);

// ============================================================================
// Bitwise Operations Benchmarks
// ============================================================================

static void BM_BigInt_LeftShift(benchmark::State& state) {
  BigInt a("123456789012345678901234567890");
  for (auto _ : state) {
    BigInt result = a << static_cast<size_t>(64);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_LeftShift);

static void BM_BigInt_RightShift(benchmark::State& state) {
  BigInt a("123456789012345678901234567890123456789012345678901234567890");
  for (auto _ : state) {
    BigInt result = a >> static_cast<size_t>(64);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_RightShift);

// ============================================================================
// HUGE Numbers Benchmarks (1000+ digits) - Real Crypto-Scale Testing!
// ============================================================================

// 1024-digit number (approximately 3400 bits - RSA-3072 scale)
static const std::string HUGE_1024 =
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "123456789012345678901234";

static const std::string HUGE_1024_B =
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "987654321098765432109876";

static void BM_BigInt_Addition_Huge(benchmark::State& state) {
  BigInt a(HUGE_1024);
  BigInt b(HUGE_1024_B);
  for (auto _ : state) {
    BigInt result = a + b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Addition_Huge);

static void BM_BigInt_Multiplication_Huge(benchmark::State& state) {
  BigInt a(HUGE_1024);
  BigInt b(HUGE_1024_B);
  for (auto _ : state) {
    BigInt result = a * b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Multiplication_Huge);

static void BM_BigInt_Division_Huge(benchmark::State& state) {
  BigInt a(HUGE_1024_B);
  BigInt b(HUGE_1024);
  for (auto _ : state) {
    BigInt result = a / b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Division_Huge);

static void BM_BigInt_GCD_Huge(benchmark::State& state) {
  BigInt a(HUGE_1024);
  BigInt b(HUGE_1024_B);
  for (auto _ : state) {
    BigInt result = gcd(a, b);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_GCD_Huge);

static void BM_BigInt_PowMod_RSA(benchmark::State& state) {
  // Simulate RSA encryption: m^e mod n with realistic sizes
  BigInt base(HUGE_1024);
  BigInt exp("65537");  // Common RSA public exponent
  BigInt mod(HUGE_1024_B);
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_RSA);

static void BM_BigInt_LeftShift_Huge(benchmark::State& state) {
  BigInt a(HUGE_1024);
  for (auto _ : state) {
    BigInt result = a << static_cast<size_t>(512);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_LeftShift_Huge);

static void BM_BigInt_RightShift_Huge(benchmark::State& state) {
  BigInt a(HUGE_1024);
  for (auto _ : state) {
    BigInt result = a >> static_cast<size_t>(512);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_RightShift_Huge);

// ============================================================================
// Crypto-Relevant Operations Benchmarks
// ============================================================================

static void BM_BigInt_IsCoprime_Small(benchmark::State& state) {
  BigInt a("12345678901234567890");
  BigInt b("98765432109876543211");  // Likely coprime
  for (auto _ : state) {
    bool result = isCoprime(a, b);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_IsCoprime_Small);

static void BM_BigInt_IsCoprime_Large(benchmark::State& state) {
  BigInt a("123456789012345678901234567890123456789012345678901234567890");
  BigInt b("987654321098765432109876543210987654321098765432109876543211");
  for (auto _ : state) {
    bool result = isCoprime(a, b);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_IsCoprime_Large);

static void BM_BigInt_InverseMod_Small(benchmark::State& state) {
  BigInt a("65537");   // Common RSA public exponent
  BigInt mod("3120");  // phi(3233)
  for (auto _ : state) {
    BigInt result = inversemod(a, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_InverseMod_Small);

static void BM_BigInt_InverseMod_Large(benchmark::State& state) {
  BigInt a("65537");
  // Large modulus (simulating phi(n) for RSA)
  BigInt mod("123456789012345678901234567890123456789012345678901234567890");
  for (auto _ : state) {
    BigInt result = inversemod(a, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_InverseMod_Large);

static void BM_BigInt_SymbolJacobi(benchmark::State& state) {
  BigInt a("12345678901234567890");
  BigInt b("98765432109876543211");  // Must be odd
  for (auto _ : state) {
    int8_t result = symbolJacobi(a, b);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_SymbolJacobi);

static void BM_BigInt_LCM_Large(benchmark::State& state) {
  BigInt a("123456789012345678901234567890");
  BigInt b("987654321098765432109876543210");
  for (auto _ : state) {
    BigInt result = lcm(a, b);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_LCM_Large);

static void BM_BigInt_Subtraction_Huge(benchmark::State& state) {
  BigInt a(HUGE_1024_B);  // Larger
  BigInt b(HUGE_1024);    // Smaller
  for (auto _ : state) {
    BigInt result = a - b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Subtraction_Huge);

static void BM_BigInt_Modulo_Huge(benchmark::State& state) {
  BigInt a(HUGE_1024_B);
  BigInt b(HUGE_1024);
  for (auto _ : state) {
    BigInt result = a % b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Modulo_Huge);

// ============================================================================
// Prime Generation and Testing Benchmarks
// ============================================================================

static void BM_BigInt_IsProbablePrime_Small(benchmark::State& state) {
  BigInt prime("104729");  // Known prime
  for (auto _ : state) {
    bool result = prime.isProbablePrime();
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_IsProbablePrime_Small);

static void BM_BigInt_IsProbablePrime_32bit(benchmark::State& state) {
  BigInt prime("2147483647");  // Mersenne prime M31
  for (auto _ : state) {
    bool result = prime.isProbablePrime();
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_IsProbablePrime_32bit);

static void BM_BigInt_RandomPrime_16bit(benchmark::State& state) {
  for (auto _ : state) {
    BigInt prime = BigInt::randomPrime(16);
    benchmark::DoNotOptimize(prime);
  }
}
BENCHMARK(BM_BigInt_RandomPrime_16bit);

static void BM_BigInt_NextPrime_Small(benchmark::State& state) {
  BigInt start("10000");
  for (auto _ : state) {
    BigInt prime = start.nextPrime();
    benchmark::DoNotOptimize(prime);
  }
}
BENCHMARK(BM_BigInt_NextPrime_Small);

static void BM_BigInt_Sqrt_Large(benchmark::State& state) {
  BigInt n("123456789012345678901234567890123456789012345678901234567890");
  for (auto _ : state) {
    BigInt result = sqrt(n);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Sqrt_Large);

static void BM_BigInt_RandomBits_256(benchmark::State& state) {
  for (auto _ : state) {
    BigInt r = BigInt::randomBits(256);
    benchmark::DoNotOptimize(r);
  }
}
BENCHMARK(BM_BigInt_RandomBits_256);

BENCHMARK_MAIN();
