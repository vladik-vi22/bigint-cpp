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

// 1024-digit ODD number for Montgomery algorithm testing (ends in 7)
static const std::string HUGE_1024_ODD =
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
    "987654321098765432109877";

// 256-digit numbers (~850 bits) for private key simulation benchmarks
// HUGE_256_ODD ends in 9 (odd) - will use Montgomery
static const std::string HUGE_256_ODD =
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567890"
    "12345678901234567890123456789012345678901234567899";

// HUGE_256_EVEN ends in 0 (even) - will use standard algorithm
static const std::string HUGE_256_EVEN =
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210"
    "98765432109876543210987654321098765432109876543210";

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

// Generate ~2600 digit string (~8640 bits, ~270 words) to trigger Toom-Cook 3-way
static std::string generateHugeString() {
  std::string result;
  for (int i = 0; i < 32; ++i) {
    result += "12345678901234567890123456789012345678901234567890123456789012345678901234567890";
  }
  return result;
}

static std::string generateHugeStringB() {
  std::string result;
  for (int i = 0; i < 32; ++i) {
    result += "98765432109876543210987654321098765432109876543210987654321098765432109876543210";
  }
  return result;
}

static void BM_BigInt_Multiplication_Toom3(benchmark::State& state) {
  // Numbers > 8192 bits (256 words) trigger Toom-Cook 3-way multiplication
  static const std::string HUGE_TOOM3_A = generateHugeString();
  static const std::string HUGE_TOOM3_B = generateHugeStringB();
  BigInt a(HUGE_TOOM3_A, 10);
  BigInt b(HUGE_TOOM3_B, 10);
  for (auto _ : state) {
    BigInt result = a * b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_Multiplication_Toom3);

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

static void BM_BigInt_PowMod_RSA_Standard(benchmark::State& state) {
  // RSA with EVEN modulus - forces standard square-and-multiply algorithm
  BigInt base(HUGE_1024);
  BigInt exp("65537");  // Common RSA public exponent (17 bits)
  BigInt mod(HUGE_1024_B);  // Even modulus (ends in 6)
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_RSA_Standard);

static void BM_BigInt_PowMod_RSA_Montgomery(benchmark::State& state) {
  // RSA with ODD modulus - uses Montgomery multiplication algorithm
  BigInt base(HUGE_1024);
  BigInt exp("65537");  // Common RSA public exponent (17 bits)
  BigInt mod(HUGE_1024_ODD);  // Odd modulus (ends in 7)
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_RSA_Montgomery);

// RSA private key simulation - large exponent where Montgomery should shine
// Uses 256-digit (~850 bit) numbers - large enough to trigger Montgomery
static void BM_BigInt_PowMod_PrivateKey_Standard(benchmark::State& state) {
  BigInt base(HUGE_256_ODD);
  BigInt exp(HUGE_256_ODD);  // Large exponent (~850 bits) - simulates private key
  BigInt mod(HUGE_256_EVEN);  // Even modulus - forces standard algorithm
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_PrivateKey_Standard);

static void BM_BigInt_PowMod_PrivateKey_Montgomery(benchmark::State& state) {
  BigInt base(HUGE_256_EVEN);
  BigInt exp(HUGE_256_ODD);  // Large exponent (~850 bits) - simulates private key
  BigInt mod(HUGE_256_ODD);  // Odd modulus - uses Montgomery
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_PrivateKey_Montgomery);

// Barrett reduction benchmarks
// Barrett: any modulus >= 4 words (128 bits), exponent >= 64 bits, not Montgomery
static void BM_BigInt_PowMod_Barrett_EvenMod(benchmark::State& state) {
  // 128-bit even modulus, 64-bit exponent -> Barrett
  BigInt base("12345678901234567890123456789012345678901234567890", 10);
  BigInt exp("18446744073709551617", 10);  // 2^64 + 1 (65 bits)
  BigInt mod("340282366920938463463374607431768211456", 10);  // 2^128 (even)
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_Barrett_EvenMod);

static void BM_BigInt_PowMod_Barrett_OddMod(benchmark::State& state) {
  // 160-bit odd modulus (< 256 bits, so not Montgomery), 65-bit exponent -> Barrett
  BigInt base("12345678901234567890123456789012345678901234567890", 10);
  BigInt exp("18446744073709551617", 10);  // 2^64 + 1 (65 bits)
  BigInt mod("1461501637330902918203684832716283019655932542983", 10);  // ~160 bits, odd
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_Barrett_OddMod);

// Compare all three algorithms with similar workload
static void BM_BigInt_PowMod_Compare_Standard(benchmark::State& state) {
  // Small modulus (< 128 bits) -> Standard
  BigInt base("12345678901234567890", 10);
  BigInt exp("1000000007", 10);  // ~30 bits
  BigInt mod("1000000007", 10);  // ~30 bits
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_Compare_Standard);

static void BM_BigInt_PowMod_Compare_Barrett(benchmark::State& state) {
  // 128-bit even modulus, 65-bit exponent -> Barrett
  BigInt base("12345678901234567890", 10);
  BigInt exp("18446744073709551617", 10);  // 65 bits
  BigInt mod("340282366920938463463374607431768211456", 10);  // 2^128
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_Compare_Barrett);

static void BM_BigInt_PowMod_Compare_Montgomery(benchmark::State& state) {
  // 256-bit odd modulus, 17-bit exponent -> Montgomery
  BigInt base("12345678901234567890", 10);
  BigInt exp("65537", 10);  // 17 bits
  BigInt mod(
      "115792089237316195423570985008687907853269984665640564039457584007913129639935",
      10);  // ~256 bits, odd
  for (auto _ : state) {
    BigInt result = powmod(base, exp, mod);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_BigInt_PowMod_Compare_Montgomery);

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
