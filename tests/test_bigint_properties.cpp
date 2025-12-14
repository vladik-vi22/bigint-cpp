/// @file test_bigint_properties.cpp
/// @brief Property-based tests for BigInt mathematical correctness.
/// @details Tests mathematical identities that must hold for all inputs.
///          Uses random values to find edge cases that unit tests miss.

#include <bigint/BigInt.hpp>

#include <gtest/gtest.h>

using namespace bigint;

/// @brief Number of random iterations per property test.
constexpr size_t kIterations = 500;

/// @brief Bit sizes to test (covers word boundaries and Karatsuba threshold).
constexpr size_t kSmallBits = 32;
constexpr size_t kMediumBits = 128;
constexpr size_t kLargeBits = 512;

class BigIntPropertyTest : public ::testing::Test {
 protected:
  /// @brief Generate random BigInt with random sign.
  static BigInt randomSigned(size_t bits) {
    BigInt value = BigInt::randomBits(bits);
    // Randomly make negative (50% chance)
    if (BigInt::randomBits(1) == BigInt(1)) {
      value = -value;
    }
    return value;
  }

  /// @brief Generate random non-zero BigInt.
  static BigInt randomNonZero(size_t bits) {
    BigInt value = BigInt::randomBits(bits);
    return value ? value : BigInt(1);
  }
};

// ============================================================================
// Addition Properties
// ============================================================================

TEST_F(BigIntPropertyTest, AdditionSubtractionIdentity) {
  // Property: (a + b) - b == a
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kMediumBits);
    BigInt b = randomSigned(kMediumBits);
    EXPECT_EQ(a + b - b, a) << "Failed for a=" << a << ", b=" << b;
  }
}

TEST_F(BigIntPropertyTest, AdditionCommutativity) {
  // Property: a + b == b + a
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kMediumBits);
    BigInt b = randomSigned(kMediumBits);
    EXPECT_EQ(a + b, b + a) << "Failed for a=" << a << ", b=" << b;
  }
}

TEST_F(BigIntPropertyTest, AdditionAssociativity) {
  // Property: (a + b) + c == a + (b + c)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kMediumBits);
    BigInt b = randomSigned(kMediumBits);
    BigInt c = randomSigned(kMediumBits);
    EXPECT_EQ((a + b) + c, a + (b + c));
  }
}

TEST_F(BigIntPropertyTest, AdditionIdentityElement) {
  // Property: a + 0 == a
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kLargeBits);
    EXPECT_EQ(a + BigInt(0), a);
    EXPECT_EQ(BigInt(0) + a, a);
  }
}

TEST_F(BigIntPropertyTest, AdditiveInverse) {
  // Property: a + (-a) == 0
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kMediumBits);
    EXPECT_EQ(a + (-a), BigInt(0));
  }
}

// ============================================================================
// Multiplication Properties
// ============================================================================

TEST_F(BigIntPropertyTest, MultiplicationDivisionIdentity) {
  // Property: (a * b) / b == a (for b != 0)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kMediumBits);
    BigInt b = randomNonZero(kSmallBits);
    if (b.isNegative()) b = -b;  // Use positive divisor for simpler test
    EXPECT_EQ((a * b) / b, a) << "Failed for a=" << a << ", b=" << b;
  }
}

TEST_F(BigIntPropertyTest, MultiplicationCommutativity) {
  // Property: a * b == b * a
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kMediumBits);
    BigInt b = randomSigned(kMediumBits);
    EXPECT_EQ(a * b, b * a);
  }
}

TEST_F(BigIntPropertyTest, MultiplicationAssociativity) {
  // Property: (a * b) * c == a * (b * c)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kSmallBits);
    BigInt b = randomSigned(kSmallBits);
    BigInt c = randomSigned(kSmallBits);
    EXPECT_EQ((a * b) * c, a * (b * c));
  }
}

TEST_F(BigIntPropertyTest, MultiplicationIdentityElement) {
  // Property: a * 1 == a
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kLargeBits);
    EXPECT_EQ(a * BigInt(1), a);
    EXPECT_EQ(BigInt(1) * a, a);
  }
}

TEST_F(BigIntPropertyTest, MultiplicationByZero) {
  // Property: a * 0 == 0
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kLargeBits);
    EXPECT_EQ(a * BigInt(0), BigInt(0));
  }
}

TEST_F(BigIntPropertyTest, DistributiveLaw) {
  // Property: a * (b + c) == a*b + a*c
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kSmallBits);
    BigInt b = randomSigned(kSmallBits);
    BigInt c = randomSigned(kSmallBits);
    EXPECT_EQ(a * (b + c), a * b + a * c);
  }
}

// ============================================================================
// Division Properties
// ============================================================================

TEST_F(BigIntPropertyTest, DivisionModuloRelation) {
  // Property: a == (a / b) * b + (a % b)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kLargeBits);
    BigInt b = randomNonZero(kMediumBits);
    EXPECT_EQ((a / b) * b + (a % b), a) << "Failed for a=" << a << ", b=" << b;
  }
}

TEST_F(BigIntPropertyTest, ModuloRange) {
  // Property: 0 <= |a % b| < |b|
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kLargeBits);
    BigInt b = randomNonZero(kMediumBits);
    BigInt remainder = a % b;
    EXPECT_LT(abs(remainder), abs(b));
  }
}

// ============================================================================
// Power and Modular Arithmetic
// ============================================================================

TEST_F(BigIntPropertyTest, PowerOfZero) {
  // Property: a^0 == 1 (for a != 0)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomNonZero(kSmallBits);
    EXPECT_EQ(pow(a, BigInt(0)), BigInt(1));
  }
}

TEST_F(BigIntPropertyTest, PowerOfOne) {
  // Property: a^1 == a
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kMediumBits);
    EXPECT_EQ(pow(a, BigInt(1)), a);
  }
}

TEST_F(BigIntPropertyTest, PowerModCorrectness) {
  // Property: powmod(a, b, m) == pow(a, b) % m (for small exponents)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt base = BigInt::randomBits(kSmallBits);
    BigInt exp = BigInt::randomBits(4);  // Small exponent to avoid huge numbers
    BigInt mod = randomNonZero(kMediumBits);
    if (mod.isNegative()) mod = -mod;
    
    BigInt expected = pow(base, exp) % mod;
    BigInt actual = powmod(base, exp, mod);
    EXPECT_EQ(actual, expected);
  }
}

// ============================================================================
// GCD and LCM Properties
// ============================================================================

TEST_F(BigIntPropertyTest, GcdDividesBoth) {
  // Property: gcd(a, b) divides both a and b
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = BigInt::randomBits(kMediumBits);
    BigInt b = BigInt::randomBits(kMediumBits);
    BigInt g = gcd(a, b);
    if (g != BigInt(0)) {
      EXPECT_EQ(a % g, BigInt(0)) << "gcd doesn't divide a";
      EXPECT_EQ(b % g, BigInt(0)) << "gcd doesn't divide b";
    }
  }
}

TEST_F(BigIntPropertyTest, GcdCommutativity) {
  // Property: gcd(a, b) == gcd(b, a)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = BigInt::randomBits(kMediumBits);
    BigInt b = BigInt::randomBits(kMediumBits);
    EXPECT_EQ(gcd(a, b), gcd(b, a));
  }
}

TEST_F(BigIntPropertyTest, LcmGcdRelation) {
  // Property: lcm(a, b) * gcd(a, b) == |a * b|
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = BigInt::randomBits(kSmallBits);
    BigInt b = BigInt::randomBits(kSmallBits);
    if (a && b) {  // Skip if either is zero
      EXPECT_EQ(lcm(a, b) * gcd(a, b), abs(a * b));
    }
  }
}

// ============================================================================
// Bitwise Properties
// ============================================================================

TEST_F(BigIntPropertyTest, LeftRightShiftInverse) {
  // Property: (a << n) >> n == a (for positive a)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = BigInt::randomBits(kMediumBits);
    size_t shift = static_cast<size_t>(BigInt::randomBits(6));  // 0-63
    EXPECT_EQ((a << shift) >> shift, a);
  }
}

TEST_F(BigIntPropertyTest, ShiftMultiplicationEquivalence) {
  // Property: a << n == a * 2^n
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = BigInt::randomBits(kMediumBits);
    size_t shift = static_cast<size_t>(BigInt::randomBits(5));  // 0-31
    BigInt power_of_two = BigInt(1) << shift;
    EXPECT_EQ(a << shift, a * power_of_two);
  }
}

TEST_F(BigIntPropertyTest, BitwiseAndSelf) {
  // Property: a & a == a
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = BigInt::randomBits(kMediumBits);
    EXPECT_EQ(a & a, a);
  }
}

TEST_F(BigIntPropertyTest, BitwiseOrSelf) {
  // Property: a | a == a
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = BigInt::randomBits(kMediumBits);
    EXPECT_EQ(a | a, a);
  }
}

TEST_F(BigIntPropertyTest, BitwiseXorSelf) {
  // Property: a ^ a == 0
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = BigInt::randomBits(kMediumBits);
    EXPECT_EQ(a ^ a, BigInt(0));
  }
}

// ============================================================================
// Comparison Properties
// ============================================================================

TEST_F(BigIntPropertyTest, ComparisonReflexivity) {
  // Property: a == a
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kLargeBits);
    EXPECT_EQ(a, a);
    EXPECT_FALSE(a < a);
    EXPECT_FALSE(a > a);
  }
}

TEST_F(BigIntPropertyTest, ComparisonAntisymmetry) {
  // Property: if a < b then !(b < a)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kMediumBits);
    BigInt b = randomSigned(kMediumBits);
    if (a < b) {
      EXPECT_FALSE(b < a);
      EXPECT_TRUE(b > a);
    }
  }
}

// ============================================================================
// String Round-Trip Properties
// ============================================================================

TEST_F(BigIntPropertyTest, DecimalStringRoundTrip) {
  // Property: BigInt(a.toStdString(10), 10) == a
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = randomSigned(kLargeBits);
    std::string str = a.toStdString(10);
    BigInt reconstructed(str, 10);
    EXPECT_EQ(reconstructed, a) << "Round-trip failed for: " << str;
  }
}

TEST_F(BigIntPropertyTest, HexStringRoundTrip) {
  // Property: BigInt(a.toStdString(16), 16) == a (for positive)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = BigInt::randomBits(kLargeBits);
    std::string str = a.toStdString(16);
    BigInt reconstructed(str, 16);
    EXPECT_EQ(reconstructed, a) << "Round-trip failed for: " << str;
  }
}

TEST_F(BigIntPropertyTest, BinaryStringRoundTrip) {
  // Property: BigInt(a.toStdString(2), 2) == a (for positive)
  for (size_t i = 0; i < kIterations; ++i) {
    BigInt a = BigInt::randomBits(kMediumBits);
    std::string str = a.toStdString(2);
    BigInt reconstructed(str, 2);
    EXPECT_EQ(reconstructed, a) << "Round-trip failed for: " << str;
  }
}

// ============================================================================
// Boundary Value Tests
// ============================================================================

TEST_F(BigIntPropertyTest, WordBoundaryOperations) {
  // Test operations around 32-bit and 64-bit boundaries
  std::vector<BigInt> boundaries = {
      BigInt(UINT32_MAX) - BigInt(1),
      BigInt(UINT32_MAX),
      BigInt(UINT32_MAX) + BigInt(1),
      BigInt(static_cast<uint64_t>(UINT32_MAX) + 1),  // 2^32
      BigInt(UINT64_MAX) - BigInt(1),
      BigInt(UINT64_MAX),
  };

  for (const auto& a : boundaries) {
    for (const auto& b : boundaries) {
      // These should not crash or produce wrong results
      BigInt sum = a + b;
      BigInt diff = a - b;
      BigInt prod = a * b;
      if (b != BigInt(0)) {
        BigInt quot = a / b;
        BigInt rem = a % b;
        EXPECT_EQ(quot * b + rem, a);
      }
    }
  }
}

TEST_F(BigIntPropertyTest, KaratsubaThresholdOperations) {
  // Test multiplication around Karatsuba threshold (~32 digits)
  for (size_t bits = 256; bits <= 1024; bits += 128) {
    BigInt a = BigInt::randomBits(bits);
    BigInt b = BigInt::randomBits(bits);
    BigInt product = a * b;
    
    // Verify via division
    if (b != BigInt(0)) {
      EXPECT_EQ(product / b, a);
    }
  }
}
