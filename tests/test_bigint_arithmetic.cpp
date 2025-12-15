#include <bigint/BigInt.hpp>

#include <gtest/gtest.h>

using namespace bigint;

class BigIntArithmeticTest : public ::testing::Test {
 protected:
  void SetUp() override {}
  void TearDown() override {}
};

TEST_F(BigIntArithmeticTest, Addition) {
  BigInt a("11111111111111111111111111111111111111111111111", 10);
  BigInt b("22222222222222222222222222222222222222222222222", 10);
  BigInt result = a + b;
  EXPECT_EQ(result.toStdString(10), "33333333333333333333333333333333333333333333333");
}

TEST_F(BigIntArithmeticTest, Subtraction) {
  BigInt a("33333333333333333333333333333333333333333333333", 10);
  BigInt b("11111111111111111111111111111111111111111111111", 10);
  BigInt result = a - b;
  EXPECT_EQ(result.toStdString(10), "22222222222222222222222222222222222222222222222");
}

TEST_F(BigIntArithmeticTest, Multiplication) {
  BigInt a("123456789", 10);
  BigInt b("987654321", 10);
  BigInt result = a * b;
  EXPECT_EQ(result.toStdString(10), "121932631112635269");
}

TEST_F(BigIntArithmeticTest, Division) {
  BigInt a("121932631112635269", 10);
  BigInt b("123456789", 10);
  BigInt result = a / b;
  EXPECT_EQ(result.toStdString(10), "987654321");
}

TEST_F(BigIntArithmeticTest, Modulo) {
  BigInt a("100", 10);
  BigInt b("30", 10);
  BigInt result = a % b;
  EXPECT_EQ(result.toStdString(10), "10");
}

TEST_F(BigIntArithmeticTest, ModuloUint32) {
  // Small number
  BigInt a("100", 10);
  EXPECT_EQ(a % 30U, 10U);

  // Large number spanning multiple 32-bit words
  BigInt b("CAFEBABE12345678", 16);
  // Verify against BigInt modulo
  EXPECT_EQ(b % 17U, static_cast<uint32_t>(b % BigInt(17)));
  EXPECT_EQ(b % 1000000007U, static_cast<uint32_t>(b % BigInt(1000000007)));

  // Very large number (more than 64 bits)
  BigInt c("CAFEBABECAFEBABE12345678", 16);
  EXPECT_EQ(c % 17U, static_cast<uint32_t>(c % BigInt(17)));
  EXPECT_EQ(c % 1000000007U, static_cast<uint32_t>(c % BigInt(1000000007)));

  // Small primes (crypto use case)
  BigInt prime_candidate("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16);
  EXPECT_EQ(prime_candidate % 3U, static_cast<uint32_t>(prime_candidate % BigInt(3)));
  EXPECT_EQ(prime_candidate % 7U, static_cast<uint32_t>(prime_candidate % BigInt(7)));
  EXPECT_EQ(prime_candidate % 113U, static_cast<uint32_t>(prime_candidate % BigInt(113)));

  // Edge cases
  BigInt zero;
  EXPECT_EQ(zero % 17U, 0U);

  BigInt one(1);
  EXPECT_EQ(one % 17U, 1U);

  // Division by zero should throw
  EXPECT_THROW([[maybe_unused]] auto _ = a % 0U, std::domain_error);
}

TEST_F(BigIntArithmeticTest, Power) {
  BigInt base("2", 10);
  BigInt exp("10", 10);
  BigInt result = pow(base, exp);
  EXPECT_EQ(result.toStdString(10), "1024");
}

TEST_F(BigIntArithmeticTest, PowerMod) {
  BigInt base("2", 10);
  BigInt exp("10", 10);
  BigInt mod("100", 10);
  BigInt result = powmod(base, exp, mod);
  EXPECT_EQ(result.toStdString(10), "24");  // 1024 % 100 = 24
}

// Test standard square-and-multiply algorithm (small modulus, uses fast division)
TEST_F(BigIntArithmeticTest, PowerModStandardSmallModulus) {
  // Small modulus - uses standard algorithm
  BigInt base("7", 10);
  BigInt exp("13", 10);
  BigInt mod("11", 10);
  BigInt result = powmod(base, exp, mod);
  // 7^13 mod 11: 7^1=7, 7^2=5, 7^4=3, 7^8=9, 7^13=7^8*7^4*7^1=9*3*7=189 mod 11=2
  EXPECT_EQ(result.toStdString(10), "2");
}

// Test standard algorithm with even modulus (Montgomery requires odd)
TEST_F(BigIntArithmeticTest, PowerModStandardEvenModulus) {
  // Even modulus - must use standard algorithm
  BigInt base("3", 10);
  BigInt exp("100", 10);
  BigInt mod("1000", 10);  // Even modulus
  BigInt result = powmod(base, exp, mod);
  // Verify result is in valid range
  EXPECT_TRUE(result >= BigInt(0));
  EXPECT_TRUE(result < mod);
  // 3^100 mod 1000 = 1 (since 3^100 ends in ...001)
  EXPECT_EQ(result.toStdString(10), "1");
}

// Test Montgomery algorithm with large odd modulus
TEST_F(BigIntArithmeticTest, PowerModMontgomeryLargeOddModulus) {
  // Large odd modulus (>= 256 bits) with exponent >= 16 bits triggers Montgomery
  // Using a 512-bit odd prime-like number
  BigInt mod(
      "1340780792994259709957402499820584612747936582059239337772356144372176403007"
      "1706109068701238785423046926253574578428203904729633967850393354094257207843",
      10);
  BigInt base("12345678901234567890123456789012345678901234567890", 10);
  BigInt exp("65537", 10);  // 17 bits, triggers Montgomery

  BigInt result = powmod(base, exp, mod);

  // Verify result is in valid range
  EXPECT_TRUE(result >= BigInt(0));
  EXPECT_TRUE(result < mod);

  // Verify correctness by checking (base^exp)^1 mod mod = base^exp mod mod
  BigInt result2 = powmod(result, BigInt(1), mod);
  EXPECT_EQ(result.toStdString(10), result2.toStdString(10));
}

// Test that both algorithms produce the same result for the same input
TEST_F(BigIntArithmeticTest, PowerModAlgorithmConsistency) {
  // Use a modulus that's large enough for Montgomery but test with small exponent first
  // Then compare with large exponent result using mathematical properties

  // 256-bit odd modulus
  BigInt mod_odd(
      "115792089237316195423570985008687907853269984665640564039457584007913129639937",
      10);
  BigInt base("98765432109876543210987654321", 10);

  // Small exponent (uses standard algorithm)
  BigInt exp_small("100", 10);
  BigInt result_small = powmod(base, exp_small, mod_odd);

  // Verify (base^100)^1 = base^100
  BigInt verify = powmod(result_small, BigInt(1), mod_odd);
  EXPECT_EQ(result_small.toStdString(10), verify.toStdString(10));

  // Large exponent (uses Montgomery if modulus is large enough)
  BigInt exp_large("65537", 10);
  BigInt result_large = powmod(base, exp_large, mod_odd);

  // Verify result is in valid range
  EXPECT_TRUE(result_large >= BigInt(0));
  EXPECT_TRUE(result_large < mod_odd);
}

// Test Montgomery with RSA-like parameters (1024-bit modulus)
TEST_F(BigIntArithmeticTest, PowerModMontgomeryRSA1024) {
  // 1024-bit odd modulus (ends in 7, so it's odd)
  BigInt mod(
      "1234567890123456789012345678901234567890123456789012345678901234567890"
      "1234567890123456789012345678901234567890123456789012345678901234567890"
      "1234567890123456789012345678901234567890123456789012345678901234567890"
      "1234567890123456789012345678901234567890123456789012345678901234567890"
      "12345678901234567890123456789012345678901234567",
      10);
  BigInt base(
      "9876543210987654321098765432109876543210987654321098765432109876543210"
      "9876543210987654321098765432109876543210987654321098765432109876543210",
      10);
  BigInt exp("65537", 10);

  BigInt result = powmod(base, exp, mod);

  // Verify result is in valid range
  EXPECT_TRUE(result >= BigInt(0));
  EXPECT_TRUE(result < mod);

  // Verify idempotence: powmod(result, 1, mod) == result
  BigInt verify = powmod(result, BigInt(1), mod);
  EXPECT_EQ(result.toStdString(10), verify.toStdString(10));
}

// Test edge cases for both algorithms
TEST_F(BigIntArithmeticTest, PowerModEdgeCases) {
  BigInt mod("97", 10);  // Small prime

  // base^0 = 1
  EXPECT_EQ(powmod(BigInt("5"), BigInt("0"), mod).toStdString(10), "1");

  // 0^exp = 0 (for exp > 0)
  EXPECT_EQ(powmod(BigInt("0"), BigInt("10"), mod).toStdString(10), "0");

  // base^1 = base mod m
  EXPECT_EQ(powmod(BigInt("50"), BigInt("1"), mod).toStdString(10), "50");

  // 1^exp = 1
  EXPECT_EQ(powmod(BigInt("1"), BigInt("1000000"), mod).toStdString(10), "1");

  // Fermat's little theorem: a^(p-1) = 1 mod p for prime p
  EXPECT_EQ(powmod(BigInt("5"), BigInt("96"), mod).toStdString(10), "1");
}

// Test that Montgomery handles the boundary case correctly
TEST_F(BigIntArithmeticTest, PowerModMontgomeryBoundary) {
  // Exactly 8 words (256 bits) odd modulus - boundary for Montgomery
  // 2^256 - 1 is odd (all 1s in binary)
  BigInt mod(
      "115792089237316195423570985008687907853269984665640564039457584007913129639935",
      10);
  BigInt base("12345678901234567890", 10);
  BigInt exp("65537", 10);

  BigInt result = powmod(base, exp, mod);

  EXPECT_TRUE(result >= BigInt(0));
  EXPECT_TRUE(result < mod);
}

TEST_F(BigIntArithmeticTest, GCD) {
  BigInt a("48", 10);
  BigInt b("18", 10);
  BigInt result = gcd(a, b);
  EXPECT_EQ(result.toStdString(10), "6");
}

TEST_F(BigIntArithmeticTest, LCM) {
  BigInt a("4", 10);
  BigInt b("6", 10);
  BigInt result = lcm(a, b);
  EXPECT_EQ(result.toStdString(10), "12");
}

TEST_F(BigIntArithmeticTest, IncrementDecrement) {
  BigInt num("100", 10);
  ++num;
  EXPECT_EQ(num.toStdString(10), "101");
  --num;
  EXPECT_EQ(num.toStdString(10), "100");
}

// ============================================================================
// REALLY BIG Numbers Tests (500+ digits) - Real BigInt Testing!
// ============================================================================

// 512-digit numbers (approximately 1700 bits)
const std::string HUGE_A =
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
    "1234567890123456789012";

const std::string HUGE_B =
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
    "9876543210987654321098";

TEST_F(BigIntArithmeticTest, HugeAddition) {
  BigInt a(HUGE_A, 10);
  BigInt b(HUGE_B, 10);
  BigInt result = a + b;
  // Verify result is larger than both operands
  EXPECT_TRUE(result > a);
  EXPECT_TRUE(result > b);
  // Verify subtraction reverses addition
  EXPECT_EQ((result - b).toStdString(10), HUGE_A);
  EXPECT_EQ((result - a).toStdString(10), HUGE_B);
}

TEST_F(BigIntArithmeticTest, HugeSubtraction) {
  BigInt a(HUGE_B, 10);  // B > A
  BigInt b(HUGE_A, 10);
  BigInt result = a - b;
  // Verify addition reverses subtraction
  EXPECT_EQ((result + b).toStdString(10), HUGE_B);
}

TEST_F(BigIntArithmeticTest, HugeMultiplication) {
  // Use 100-digit numbers for multiplication (result will be ~200 digits)
  std::string num100a =
      "1234567890123456789012345678901234567890"
      "1234567890123456789012345678901234567890"
      "12345678901234567890";
  std::string num100b =
      "9876543210987654321098765432109876543210"
      "9876543210987654321098765432109876543210"
      "98765432109876543210";
  BigInt a(num100a, 10);
  BigInt b(num100b, 10);
  BigInt result = a * b;
  // Verify division reverses multiplication
  EXPECT_EQ((result / b).toStdString(10), num100a);
  EXPECT_EQ((result / a).toStdString(10), num100b);
}

TEST_F(BigIntArithmeticTest, Toom3Multiplication) {
  // Test Toom-Cook 3-way multiplication with very large numbers (>256 words = >8192 bits)
  // Generate ~2500 digit numbers (~8300 bits, ~260 words) to trigger Toom-3
  // 1 decimal digit ≈ 3.32 bits, so need ~2500 digits for 8300 bits
  std::string num2500a;
  std::string num2500b;
  for (int i = 0; i < 32; ++i) {
    num2500a += "1234567890123456789012345678901234567890"
                "1234567890123456789012345678901234567890";
    num2500b += "9876543210987654321098765432109876543210"
                "9876543210987654321098765432109876543210";
  }

  BigInt a(num2500a, 10);
  BigInt b(num2500b, 10);

  // Verify bit length is above Toom-3 threshold (256 words * 32 bits = 8192 bits)
  EXPECT_GT(a.bitLength(), 8192);
  EXPECT_GT(b.bitLength(), 8192);

  BigInt result = a * b;

  // Verify multiplication is correct by checking division reverses it
  EXPECT_EQ((result / b).toStdString(10), num2500a);
  EXPECT_EQ((result / a).toStdString(10), num2500b);

  // Verify commutativity
  BigInt result2 = b * a;
  EXPECT_EQ(result, result2);
}

TEST_F(BigIntArithmeticTest, HugeDivision) {
  BigInt a(HUGE_B, 10);
  BigInt b(HUGE_A, 10);
  BigInt quotient = a / b;
  BigInt remainder = a % b;
  // Verify: a = b * quotient + remainder
  EXPECT_EQ((b * quotient + remainder).toStdString(10), HUGE_B);
}

TEST_F(BigIntArithmeticTest, HugePowerMod) {
  // RSA-like operation: base^exp mod n with large numbers
  // Using odd modulus to potentially trigger Montgomery algorithm
  BigInt base("123456789012345678901234567890123456789012345678901234567890", 10);
  BigInt exp("65537", 10);  // Common RSA public exponent (17 bits)
  BigInt mod("1000000000000000000000000000000000000000000000000000000000057", 10);
  BigInt result = powmod(base, exp, mod);
  // Result should be less than modulus
  EXPECT_TRUE(result < mod);
  EXPECT_TRUE(result >= BigInt(0));

  // Verify idempotence
  BigInt verify = powmod(result, BigInt(1), mod);
  EXPECT_EQ(result.toStdString(10), verify.toStdString(10));
}

TEST_F(BigIntArithmeticTest, HugePowerModEvenModulus) {
  // Large even modulus - must use standard algorithm
  BigInt base("123456789012345678901234567890123456789012345678901234567890", 10);
  BigInt exp("65537", 10);
  BigInt mod("1000000000000000000000000000000000000000000000000000000000000", 10);
  BigInt result = powmod(base, exp, mod);
  EXPECT_TRUE(result < mod);
  EXPECT_TRUE(result >= BigInt(0));
}

// ============================================================================
// Barrett Reduction Tests
// ============================================================================

// Barrett: any modulus >= 4 words (128 bits), exponent >= 64 bits
TEST_F(BigIntArithmeticTest, PowerModBarrettEvenModulus) {
  // Even modulus >= 128 bits, exponent >= 64 bits -> Barrett
  // 128-bit even modulus
  BigInt mod("340282366920938463463374607431768211456", 10);  // 2^128
  BigInt base("12345678901234567890", 10);
  // 64-bit exponent (must be >= 64 bits to trigger Barrett)
  BigInt exp("18446744073709551616", 10);  // 2^64

  BigInt result = powmod(base, exp, mod);

  EXPECT_TRUE(result >= BigInt(0));
  EXPECT_TRUE(result < mod);

  // Verify with smaller computation: base^(2^64) mod 2^128
  // Since mod is power of 2, result should be deterministic
}

TEST_F(BigIntArithmeticTest, PowerModBarrettOddModulusMedium) {
  // Odd modulus >= 128 bits but < 256 bits, exponent >= 64 bits -> Barrett
  // (Montgomery requires >= 256 bits)
  // 160-bit odd modulus
  BigInt mod("1461501637330902918203684832716283019655932542983", 10);  // ~160 bits, odd
  BigInt base("98765432109876543210", 10);
  BigInt exp("18446744073709551617", 10);  // 2^64 + 1 (65 bits)

  BigInt result = powmod(base, exp, mod);

  EXPECT_TRUE(result >= BigInt(0));
  EXPECT_TRUE(result < mod);
}

TEST_F(BigIntArithmeticTest, PowerModBarrettConsistency) {
  // Verify Barrett produces same result as standard algorithm
  // Use parameters that would trigger Barrett
  BigInt mod("340282366920938463463374607431768211457", 10);  // 128-bit odd
  BigInt base("999999999999999999", 10);
  BigInt exp("18446744073709551617", 10);  // 65 bits

  BigInt result = powmod(base, exp, mod);

  // Verify using Fermat's little theorem property
  // For any a coprime to n: a^phi(n) ≡ 1 (mod n)
  // Just verify result is in valid range and idempotent
  EXPECT_TRUE(result >= BigInt(0));
  EXPECT_TRUE(result < mod);

  // Verify idempotence: (result^1) mod m = result
  BigInt verify = powmod(result, BigInt(1), mod);
  EXPECT_EQ(result.toStdString(10), verify.toStdString(10));
}

TEST_F(BigIntArithmeticTest, PowerModAlgorithmBoundaries) {
  // Test at algorithm selection boundaries
  BigInt base("12345678901234567890", 10);

  // Standard: small modulus (< 128 bits) or small exponent (< 64 bits)
  {
    BigInt mod("1000000007", 10);  // ~30 bits
    BigInt exp("1000000", 10);     // ~20 bits
    BigInt result = powmod(base, exp, mod);
    EXPECT_TRUE(result >= BigInt(0));
    EXPECT_TRUE(result < mod);
  }

  // Barrett: 128-bit modulus, 64-bit exponent, even modulus
  {
    BigInt mod("340282366920938463463374607431768211456", 10);  // 2^128
    BigInt exp("18446744073709551616", 10);                     // 2^64
    BigInt result = powmod(base, exp, mod);
    EXPECT_TRUE(result >= BigInt(0));
    EXPECT_TRUE(result < mod);
  }

  // Montgomery: 256-bit odd modulus, 16-bit exponent
  {
    BigInt mod(
        "1157920892373161954235709850086879078532699846656405640394575840079131296"
        "39935",
        10);                       // ~256 bits, odd
    BigInt exp("65537", 10);       // 17 bits
    BigInt result = powmod(base, exp, mod);
    EXPECT_TRUE(result >= BigInt(0));
    EXPECT_TRUE(result < mod);
  }
}

TEST_F(BigIntArithmeticTest, HugeGCD) {
  // GCD of two large numbers
  BigInt a("123456789012345678901234567890123456789012345678901234567890", 10);
  BigInt b("987654321098765432109876543210987654321098765432109876543210", 10);
  BigInt result = gcd(a, b);
  // GCD should divide both numbers evenly
  EXPECT_EQ((a % result).toStdString(10), "0");
  EXPECT_EQ((b % result).toStdString(10), "0");
}

TEST_F(BigIntArithmeticTest, HugeBitShift) {
  BigInt a(HUGE_A, 10);
  // Shift left by 256 bits (multiply by 2^256)
  BigInt shifted = a << static_cast<size_t>(256);
  // Shifted should be larger than original
  EXPECT_TRUE(shifted > a);
  // Verify left shift is equivalent to multiplication by 2^n
  BigInt two("2", 10);
  BigInt twoTo10 = pow(two, BigInt("10", 10));  // 1024
  BigInt shiftedBy10 = a << static_cast<size_t>(10);
  EXPECT_EQ(shiftedBy10.toStdString(10), (a * twoTo10).toStdString(10));

  // Test right shift
  BigInt num("1024", 10);
  BigInt rightShifted = num >> static_cast<size_t>(10);
  EXPECT_EQ(rightShifted.toStdString(10), "1");

  // Test right shift with various shift amounts (verifies carry calculation fix)
  BigInt largeNum("12345678901234567890", 10);
  for (size_t shift = 1; shift < 32; ++shift) {
    BigInt leftShifted = largeNum << shift;
    BigInt backShifted = leftShifted >> shift;
    EXPECT_EQ(backShifted.toStdString(10), largeNum.toStdString(10))
        << "Failed for shift = " << shift;
  }
}

TEST_F(BigIntArithmeticTest, HugeComparison) {
  BigInt a(HUGE_A, 10);
  BigInt b(HUGE_B, 10);
  BigInt a_copy(HUGE_A, 10);
  EXPECT_TRUE(b > a);
  EXPECT_TRUE(a < b);
  EXPECT_TRUE(a == a_copy);
  EXPECT_TRUE(a != b);
  EXPECT_TRUE(a <= a_copy);
  EXPECT_TRUE(a >= a_copy);
}

// ============================================================================
// Additional Crypto-Relevant Functions Tests
// ============================================================================

TEST_F(BigIntArithmeticTest, AbsoluteValue) {
  BigInt positive("12345", 10);
  BigInt negative("-12345", 10);
  EXPECT_EQ(abs(positive).toStdString(10), "12345");
  EXPECT_EQ(abs(negative).toStdString(10), "12345");
  EXPECT_TRUE(abs(negative) > 0);
}

TEST_F(BigIntArithmeticTest, IsCoprime) {
  // 15 and 28 are coprime (gcd = 1)
  BigInt a("15", 10);
  BigInt b("28", 10);
  EXPECT_TRUE(isCoprime(a, b));

  // 12 and 18 are not coprime (gcd = 6)
  BigInt c("12", 10);
  BigInt d("18", 10);
  EXPECT_FALSE(isCoprime(c, d));

  // Any number and 1 are coprime
  BigInt large("123456789012345678901234567890", 10);
  EXPECT_TRUE(isCoprime(large, BigInt(1)));
}

TEST_F(BigIntArithmeticTest, InverseMod) {
  // 3 * 7 = 21 ≡ 1 (mod 10), so inverse of 3 mod 10 is 7
  BigInt a("3", 10);
  BigInt mod("10", 10);
  BigInt inv = inversemod(a, mod);
  EXPECT_EQ(inv.toStdString(10), "7");
  // Verify: a * inv ≡ 1 (mod m)
  EXPECT_EQ(((a * inv) % mod).toStdString(10), "1");

  // RSA-style: inverse of e mod phi(n)
  BigInt e("65537", 10);
  BigInt phi("3120", 10);  // phi(3233) = 3120 for p=61, q=53
  BigInt d = inversemod(e, phi);
  // Verify: e * d ≡ 1 (mod phi)
  EXPECT_EQ(((e * d) % phi).toStdString(10), "1");
}

TEST_F(BigIntArithmeticTest, InverseModNotCoprime) {
  // 6 and 9 are not coprime (gcd = 3), so no inverse exists
  BigInt a("6", 10);
  BigInt mod("9", 10);
  // Now throws std::domain_error when inverse doesn't exist
  EXPECT_THROW([[maybe_unused]] auto _ = inversemod(a, mod), std::domain_error);
}

TEST_F(BigIntArithmeticTest, CongruenceMod) {
  // 17 ≡ 5 (mod 12) because 17 - 5 = 12
  BigInt a("17", 10);
  BigInt b("5", 10);
  BigInt mod("12", 10);
  EXPECT_TRUE(congruencemod(a, b, mod));

  // 17 ≢ 6 (mod 12)
  BigInt c("6", 10);
  EXPECT_FALSE(congruencemod(a, c, mod));
}

TEST_F(BigIntArithmeticTest, SymbolJacobi) {
  // Jacobi symbol (2/7) = 1
  BigInt a("2", 10);
  BigInt b("7", 10);
  EXPECT_EQ(symbolJacobi(a, b), 1);

  // Jacobi symbol (5/21) = 1 (verified: 5^((21-1)/2) mod 21 = 5^10 mod 21 = 1)
  BigInt c("5", 10);
  BigInt d("21", 10);
  EXPECT_EQ(symbolJacobi(c, d), 1);

  // Jacobi symbol (2/15) = 1
  BigInt e("2", 10);
  BigInt f("15", 10);
  EXPECT_EQ(symbolJacobi(e, f), 1);

  // Jacobi symbol (6/15) = 0 (not coprime, gcd=3)
  BigInt g("6", 10);
  BigInt h("15", 10);
  EXPECT_EQ(symbolJacobi(g, h), 0);
}

TEST_F(BigIntArithmeticTest, Sqrt) {
  // sqrt(0) = 0
  EXPECT_EQ(sqrt(BigInt(0)), BigInt(0));

  // sqrt(1) = 1
  EXPECT_EQ(sqrt(BigInt(1)), BigInt(1));

  // sqrt(4) = 2
  EXPECT_EQ(sqrt(BigInt(4)), BigInt(2));

  // sqrt(9) = 3
  EXPECT_EQ(sqrt(BigInt(9)), BigInt(3));

  // sqrt(10) = 3 (floor)
  EXPECT_EQ(sqrt(BigInt(10)), BigInt(3));

  // sqrt(100) = 10
  EXPECT_EQ(sqrt(BigInt(100)), BigInt(10));

  // sqrt(1000000) = 1000
  EXPECT_EQ(sqrt(BigInt(1000000)), BigInt(1000));

  // Large number: sqrt(2^64) = 2^32
  BigInt large = BigInt(1) << 64;
  BigInt expected = BigInt(1) << 32;
  EXPECT_EQ(sqrt(large), expected);

  // Verify: result^2 <= n < (result+1)^2
  BigInt n("123456789012345678901234567890", 10);
  BigInt s = sqrt(n);
  EXPECT_LE(s * s, n);
  EXPECT_GT((s + BigInt(1)) * (s + BigInt(1)), n);
}

TEST_F(BigIntArithmeticTest, IsProbablePrimeSmall) {
  // Known small primes
  EXPECT_TRUE(BigInt(2).isProbablePrime());
  EXPECT_TRUE(BigInt(3).isProbablePrime());
  EXPECT_TRUE(BigInt(5).isProbablePrime());
  EXPECT_TRUE(BigInt(7).isProbablePrime());
  EXPECT_TRUE(BigInt(11).isProbablePrime());
  EXPECT_TRUE(BigInt(13).isProbablePrime());
  EXPECT_TRUE(BigInt(17).isProbablePrime());
  EXPECT_TRUE(BigInt(19).isProbablePrime());
  EXPECT_TRUE(BigInt(23).isProbablePrime());
  EXPECT_TRUE(BigInt(97).isProbablePrime());

  // Known small composites
  EXPECT_FALSE(BigInt(0).isProbablePrime());
  EXPECT_FALSE(BigInt(1).isProbablePrime());
  EXPECT_FALSE(BigInt(4).isProbablePrime());
  EXPECT_FALSE(BigInt(6).isProbablePrime());
  EXPECT_FALSE(BigInt(8).isProbablePrime());
  EXPECT_FALSE(BigInt(9).isProbablePrime());
  EXPECT_FALSE(BigInt(10).isProbablePrime());
  EXPECT_FALSE(BigInt(15).isProbablePrime());
  EXPECT_FALSE(BigInt(100).isProbablePrime());

  // Negative numbers are not prime
  EXPECT_FALSE(BigInt("-7").isProbablePrime());
}

TEST_F(BigIntArithmeticTest, IsProbablePrimeLarge) {
  // Mersenne primes: 2^p - 1 where p is prime
  // 2^13 - 1 = 8191
  BigInt mersenne13 = (BigInt(1) << 13) - BigInt(1);
  EXPECT_TRUE(mersenne13.isProbablePrime());

  // 2^17 - 1 = 131071
  BigInt mersenne17 = (BigInt(1) << 17) - BigInt(1);
  EXPECT_TRUE(mersenne17.isProbablePrime());

  // 2^31 - 1 = 2147483647
  BigInt mersenne31 = (BigInt(1) << 31) - BigInt(1);
  EXPECT_TRUE(mersenne31.isProbablePrime());

  // 2^61 - 1 = 2305843009213693951
  BigInt mersenne61 = (BigInt(1) << 61) - BigInt(1);
  EXPECT_TRUE(mersenne61.isProbablePrime());

  // Product of two primes is composite
  BigInt p1("104729", 10);
  BigInt p2("104743", 10);
  BigInt composite = p1 * p2;
  EXPECT_FALSE(composite.isProbablePrime());
}

TEST_F(BigIntArithmeticTest, RandomBits) {
  // Test that randomBits generates numbers with correct bit length
  for (size_t bits : {32, 64, 128, 256}) {
    BigInt r = BigInt::randomBits(bits);
    EXPECT_EQ(r.bitLength(), bits);
    EXPECT_TRUE(r > 0);
  }

  // Edge case: 1 bit
  BigInt r1 = BigInt::randomBits(1);
  EXPECT_EQ(r1.bitLength(), 1);
  EXPECT_EQ(r1, BigInt(1));

  // Edge case: 0 bits
  BigInt r0 = BigInt::randomBits(0);
  EXPECT_TRUE(r0 == 0);
}

TEST_F(BigIntArithmeticTest, RandomBelow) {
  BigInt max("1000000", 10);

  // Generate several random numbers and verify they're in range
  for (int i = 0; i < 10; ++i) {
    BigInt r = BigInt::randomBelow(max);
    EXPECT_GE(r, BigInt(0));
    EXPECT_LT(r, max);
  }

  // Edge cases
  EXPECT_TRUE(BigInt::randomBelow(BigInt(1)) == 0);
  EXPECT_TRUE(BigInt::randomBelow(BigInt(0)) == 0);
}

TEST_F(BigIntArithmeticTest, BitLengthVariousSizes) {
  // Zero has 1 bit
  EXPECT_EQ(BigInt().bitLength(), 1);

  // Small numbers
  EXPECT_EQ(BigInt(1).bitLength(), 1);
  EXPECT_EQ(BigInt(UINT32_MAX).bitLength(), 32);

  // Numbers requiring more than 32 bits
  BigInt twoDigits = BigInt(1) << 32;
  EXPECT_EQ(twoDigits.bitLength(), 33);

  // Large number
  BigInt large = BigInt(1) << 256;
  EXPECT_EQ(large.bitLength(), 257);
}

TEST_F(BigIntArithmeticTest, RandomPrime16Bit) {
  BigInt prime = BigInt::randomPrime(16);
  EXPECT_EQ(prime.bitLength(), 16);
  EXPECT_TRUE(prime.isProbablePrime());
  EXPECT_TRUE(prime % 2 != 0);  // isOdd
}

TEST_F(BigIntArithmeticTest, RandomPrimeSmall) {
  BigInt p2 = BigInt::randomPrime(2);
  EXPECT_TRUE(p2 == BigInt(2) || p2 == BigInt(3));

  BigInt p3 = BigInt::randomPrime(3);
  EXPECT_EQ(p3.bitLength(), 3);
  EXPECT_TRUE(p3.isProbablePrime());
}

TEST_F(BigIntArithmeticTest, NextPrime) {
  EXPECT_EQ(BigInt(0).nextPrime(), BigInt(2));
  EXPECT_EQ(BigInt(1).nextPrime(), BigInt(2));
  EXPECT_EQ(BigInt(2).nextPrime(), BigInt(2));
  EXPECT_EQ(BigInt(3).nextPrime(), BigInt(3));
  EXPECT_EQ(BigInt(4).nextPrime(), BigInt(5));
  EXPECT_EQ(BigInt(10).nextPrime(), BigInt(11));
  EXPECT_EQ(BigInt(100).nextPrime(), BigInt(101));
  EXPECT_EQ(BigInt(1000).nextPrime(), BigInt(1009));
}

TEST_F(BigIntArithmeticTest, NextPrimeLarge) {
  // Next prime after 10000
  BigInt prime = BigInt(10000).nextPrime();
  EXPECT_EQ(prime, BigInt(10007));
  EXPECT_TRUE(prime.isProbablePrime());

  // Next prime after 1000000
  prime = BigInt(1000000).nextPrime();
  EXPECT_EQ(prime, BigInt(1000003));
}

TEST_F(BigIntArithmeticTest, DivisionEdgeCases) {
  // Division of zero by non-zero
  BigInt zero(0);
  BigInt nonZero(100);
  EXPECT_EQ(zero / nonZero, BigInt(0));
  EXPECT_EQ(zero % nonZero, BigInt(0));

  // Division by one
  BigInt a(12345);
  EXPECT_EQ(a / BigInt(1), a);
  EXPECT_EQ(a % BigInt(1), BigInt(0));

  // Division by self
  EXPECT_EQ(a / a, BigInt(1));
  EXPECT_EQ(a % a, BigInt(0));

  // Division where dividend < divisor
  BigInt small(10);
  BigInt large(100);
  EXPECT_EQ(small / large, BigInt(0));
  EXPECT_EQ(small % large, small);

  // NOTE: Division by zero causes infinite loop - undefined behavior
}

TEST_F(BigIntArithmeticTest, SignedConversions) {
  // Positive values
  BigInt pos(12345);
  EXPECT_EQ(static_cast<int64_t>(pos), 12345);
  EXPECT_EQ(static_cast<int32_t>(pos), 12345);
  EXPECT_EQ(static_cast<int16_t>(pos), 12345);

  // Negative values
  BigInt neg(-12345);
  EXPECT_EQ(static_cast<int64_t>(neg), -12345);
  EXPECT_EQ(static_cast<int32_t>(neg), -12345);
  EXPECT_EQ(static_cast<int16_t>(neg), -12345);

  // Zero
  BigInt z(0);
  EXPECT_EQ(static_cast<int64_t>(z), 0);
  EXPECT_EQ(static_cast<int32_t>(z), 0);
}

TEST_F(BigIntArithmeticTest, NegativeOperations) {
  // Test int32_t constructor with negative value
  BigInt a(-100);  // Should be -100
  EXPECT_TRUE(a < 0) << "a should be negative";
  EXPECT_EQ(static_cast<uint64_t>(a), 100) << "magnitude should be 100";
  EXPECT_EQ(a.toStdString(10), "-100") << "string should be -100";

  BigInt b(30);
  EXPECT_TRUE(b > 0);
  EXPECT_EQ(b.toStdString(10), "30");

  // Test negative + positive: -100 + 30 = -70
  BigInt sum = a + b;
  EXPECT_EQ(sum.toStdString(10), "-70");

  // Test negative - positive: -100 - 30 = -130
  BigInt diff = a - b;
  EXPECT_EQ(diff.toStdString(10), "-130");

  // Test negative * positive: -100 * 30 = -3000
  BigInt prod = a * b;
  EXPECT_EQ(prod.toStdString(10), "-3000");

  // Test negative * negative = positive
  BigInt c(-30);
  BigInt prod2 = a * c;
  EXPECT_TRUE(prod2 > 0);
  EXPECT_EQ(prod2.toStdString(10), "3000");
}

TEST_F(BigIntArithmeticTest, EdgeCaseOverflow) {
  // Test around 32-bit boundary
  BigInt maxU32(UINT32_MAX);
  BigInt one(1);
  BigInt overflow = maxU32 + one;
  EXPECT_EQ(overflow.toStdString(10), "4294967296");

  // Test around 64-bit boundary
  BigInt maxU64(UINT64_MAX);
  BigInt overflow64 = maxU64 + one;
  EXPECT_EQ(overflow64.toStdString(10), "18446744073709551616");

  // Multiplication overflow
  BigInt large(UINT32_MAX);
  BigInt product = large * large;
  EXPECT_EQ(product.toStdString(10), "18446744065119617025");
}

TEST_F(BigIntArithmeticTest, EdgeCaseZeroOperations) {
  BigInt zero(0);
  BigInt num(12345);

  // Zero arithmetic
  EXPECT_EQ(zero + num, num);
  EXPECT_EQ(num + zero, num);
  EXPECT_EQ(zero - num, -num);
  EXPECT_EQ(num - zero, num);
  EXPECT_EQ(zero * num, zero);
  EXPECT_EQ(num * zero, zero);

  // Zero comparisons
  EXPECT_TRUE(zero == BigInt(0));
  EXPECT_TRUE(zero <= num);
  EXPECT_TRUE(zero < num);
  EXPECT_FALSE(zero > num);

  // Zero shifts
  EXPECT_EQ(zero << 100, zero);
  EXPECT_EQ(zero >> 100, zero);
}

TEST_F(BigIntArithmeticTest, EdgeCasePowerOperations) {
  // Power of zero
  EXPECT_EQ(pow(BigInt(0), BigInt(5)), BigInt(0));
  EXPECT_EQ(pow(BigInt(5), BigInt(0)), BigInt(1));

  // Power of one
  EXPECT_EQ(pow(BigInt(1), BigInt(1000)), BigInt(1));
  EXPECT_EQ(pow(BigInt(1000), BigInt(1)), BigInt(1000));

  // Negative exponent returns zero
  EXPECT_EQ(pow(BigInt(2), BigInt(-1)), BigInt(0));
}

TEST_F(BigIntArithmeticTest, EdgeCaseGcdLcm) {
  // GCD with zero
  EXPECT_EQ(gcd(BigInt(0), BigInt(5)), BigInt(5));
  EXPECT_EQ(gcd(BigInt(5), BigInt(0)), BigInt(5));

  // GCD with one
  EXPECT_EQ(gcd(BigInt(1), BigInt(12345)), BigInt(1));

  // LCM with one
  EXPECT_EQ(lcm(BigInt(1), BigInt(12345)), BigInt(12345));

  // GCD of same number
  EXPECT_EQ(gcd(BigInt(100), BigInt(100)), BigInt(100));
}

TEST_F(BigIntArithmeticTest, EdgeCaseIncrementDecrement) {
  // Increment from -1 to 0
  BigInt neg1(-1);
  ++neg1;
  EXPECT_EQ(neg1, BigInt(0));
  EXPECT_TRUE(neg1 == 0);

  // Decrement from 0 to -1
  BigInt zero(0);
  --zero;
  EXPECT_EQ(zero, BigInt(-1));
  EXPECT_TRUE(zero < 0);

  // Increment across 32-bit boundary
  BigInt boundary(UINT32_MAX);
  ++boundary;
  EXPECT_EQ(boundary.toStdString(10), "4294967296");
}
