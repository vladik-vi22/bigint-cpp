#include <bigint/BigInt.hpp>

#include <gtest/gtest.h>

using namespace bigint;

class BigIntBitwiseTest : public ::testing::Test {
 protected:
  void SetUp() override {}
  void TearDown() override {}
};

TEST_F(BigIntBitwiseTest, BitwiseAND) {
  BigInt a("15", 10);  // 1111
  BigInt b("9", 10);   // 1001
  BigInt result = a & b;
  EXPECT_EQ(result.toStdString(10), "9");  // 1001
}

TEST_F(BigIntBitwiseTest, BitwiseOR) {
  BigInt a("12", 10);  // 1100
  BigInt b("10", 10);  // 1010
  BigInt result = a | b;
  EXPECT_EQ(result.toStdString(10), "14");  // 1110
}

TEST_F(BigIntBitwiseTest, BitwiseXOR) {
  BigInt a("12", 10);  // 1100
  BigInt b("10", 10);  // 1010
  BigInt result = a ^ b;
  EXPECT_EQ(result.toStdString(10), "6");  // 0110
}

TEST_F(BigIntBitwiseTest, LeftShift) {
  BigInt num("1", 10);
  BigInt result = num << static_cast<size_t>(10);
  EXPECT_EQ(result.toStdString(10), "1024");
}

TEST_F(BigIntBitwiseTest, RightShift) {
  BigInt num("1024", 10);
  BigInt result = num >> static_cast<size_t>(10);
  EXPECT_EQ(result.toStdString(10), "1");
}

// ============================================================================
// Shift Edge Cases
// ============================================================================

TEST_F(BigIntBitwiseTest, ShiftByZero) {
  BigInt num("12345", 10);
  EXPECT_EQ((num << static_cast<size_t>(0)).toStdString(10), "12345");
  EXPECT_EQ((num >> static_cast<size_t>(0)).toStdString(10), "12345");
}

TEST_F(BigIntBitwiseTest, ShiftZero) {
  BigInt zero;
  EXPECT_TRUE((zero << static_cast<size_t>(10)) == 0);
  EXPECT_TRUE((zero >> static_cast<size_t>(10)) == 0);
}

TEST_F(BigIntBitwiseTest, RightShiftToZero) {
  BigInt num("255", 10);
  BigInt result = num >> static_cast<size_t>(100);  // Shift more than bit length
  EXPECT_TRUE(result == 0);
}

TEST_F(BigIntBitwiseTest, LargeShift) {
  BigInt num("1", 10);
  BigInt shifted = num << static_cast<size_t>(1000);  // 2^1000
  EXPECT_EQ(shifted.bitLength(), 1001);

  // Shift back
  BigInt back = shifted >> static_cast<size_t>(1000);
  EXPECT_EQ(back.toStdString(10), "1");
}

TEST_F(BigIntBitwiseTest, BitwiseWithZero) {
  BigInt num("12345", 10);
  BigInt zero;

  EXPECT_TRUE((num & zero) == 0);
  EXPECT_EQ((num | zero).toStdString(10), "12345");
  EXPECT_EQ((num ^ zero).toStdString(10), "12345");
}

TEST_F(BigIntBitwiseTest, BitwiseWithSelf) {
  BigInt num("12345", 10);

  EXPECT_EQ((num & num).toStdString(10), "12345");
  EXPECT_EQ((num | num).toStdString(10), "12345");
  EXPECT_TRUE((num ^ num) == 0);
}

// ============================================================================
// testBit Tests
// ============================================================================

TEST_F(BigIntBitwiseTest, TestBitBasic) {
  BigInt num("13", 10);  // Binary: 1101
  EXPECT_TRUE(num.testBit(0));   // bit 0 = 1
  EXPECT_FALSE(num.testBit(1));  // bit 1 = 0
  EXPECT_TRUE(num.testBit(2));   // bit 2 = 1
  EXPECT_TRUE(num.testBit(3));   // bit 3 = 1
  EXPECT_FALSE(num.testBit(4));  // bit 4 = 0 (beyond value)
}

TEST_F(BigIntBitwiseTest, TestBitPowerOfTwo) {
  BigInt num = BigInt(1) << 100;  // 2^100
  EXPECT_FALSE(num.testBit(99));
  EXPECT_TRUE(num.testBit(100));
  EXPECT_FALSE(num.testBit(101));
}

TEST_F(BigIntBitwiseTest, TestBitZero) {
  BigInt zero;
  EXPECT_FALSE(zero.testBit(0));
  EXPECT_FALSE(zero.testBit(100));
  EXPECT_FALSE(zero.testBit(1000));
}

TEST_F(BigIntBitwiseTest, TestBitLargeNumber) {
  // 2^256 - 1 (all bits set for 256 bits)
  BigInt num(
      "115792089237316195423570985008687907853269984665640564039457584007913129639935",
      10);
  // All bits 0-255 should be set
  for (size_t i = 0; i < 256; ++i) {
    EXPECT_TRUE(num.testBit(i)) << "Bit " << i << " should be set";
  }
  // Bit 256 should not be set
  EXPECT_FALSE(num.testBit(256));
}

TEST_F(BigIntBitwiseTest, TestBitWordBoundary) {
  // Test bits at word boundaries (32-bit words)
  BigInt num = BigInt(1) << 32;  // 2^32
  EXPECT_FALSE(num.testBit(31));
  EXPECT_TRUE(num.testBit(32));
  EXPECT_FALSE(num.testBit(33));

  BigInt num2 = BigInt(1) << 64;  // 2^64
  EXPECT_FALSE(num2.testBit(63));
  EXPECT_TRUE(num2.testBit(64));
  EXPECT_FALSE(num2.testBit(65));
}

TEST_F(BigIntBitwiseTest, TestBitOutOfRange) {
  BigInt num("255", 10);  // 8 bits
  // Bits beyond the number's size should return false
  EXPECT_FALSE(num.testBit(100));
  EXPECT_FALSE(num.testBit(1000));
  EXPECT_FALSE(num.testBit(10000));
}

// ============================================================================
// popCount Tests
// ============================================================================

TEST_F(BigIntBitwiseTest, PopCountBasic) {
  EXPECT_EQ(BigInt(0).popCount(), 0);
  EXPECT_EQ(BigInt(1).popCount(), 1);
  EXPECT_EQ(BigInt(5).popCount(), 2);   // 101
  EXPECT_EQ(BigInt(7).popCount(), 3);   // 111
  EXPECT_EQ(BigInt(8).popCount(), 1);   // 1000
  EXPECT_EQ(BigInt(255).popCount(), 8); // 11111111
}

TEST_F(BigIntBitwiseTest, PopCountLarge) {
  // 2^256 - 1 (all 256 bits set)
  BigInt allOnes(
      "115792089237316195423570985008687907853269984665640564039457584007913129639935",
      10);
  EXPECT_EQ(allOnes.popCount(), 256);

  // 2^100 (single bit set)
  BigInt powerOfTwo = BigInt(1) << 100;
  EXPECT_EQ(powerOfTwo.popCount(), 1);
}

TEST_F(BigIntBitwiseTest, PopCountWordBoundary) {
  // 2^32 - 1 (32 bits set)
  BigInt num(0xFFFFFFFFULL);
  EXPECT_EQ(num.popCount(), 32);

  // 2^64 - 1 (64 bits set)
  BigInt num64("18446744073709551615", 10);
  EXPECT_EQ(num64.popCount(), 64);
}

// ============================================================================
// setBit Tests
// ============================================================================

TEST_F(BigIntBitwiseTest, SetBitBasic) {
  BigInt num(0);
  num.setBit(0);
  EXPECT_EQ(num, BigInt(1));

  num.setBit(2);
  EXPECT_EQ(num, BigInt(5));  // 101

  num.setBit(1);
  EXPECT_EQ(num, BigInt(7));  // 111
}

TEST_F(BigIntBitwiseTest, SetBitAlreadySet) {
  BigInt num(7);  // 111
  num.setBit(0);
  EXPECT_EQ(num, BigInt(7));  // Still 111
}

TEST_F(BigIntBitwiseTest, SetBitExtends) {
  BigInt num(1);
  num.setBit(100);
  EXPECT_TRUE(num.testBit(0));
  EXPECT_TRUE(num.testBit(100));
  EXPECT_EQ(num.bitLength(), 101);
}

TEST_F(BigIntBitwiseTest, SetBitChaining) {
  BigInt num(0);
  num.setBit(0).setBit(1).setBit(2);
  EXPECT_EQ(num, BigInt(7));
}

// ============================================================================
// clearBit Tests
// ============================================================================

TEST_F(BigIntBitwiseTest, ClearBitBasic) {
  BigInt num(7);  // 111
  num.clearBit(1);
  EXPECT_EQ(num, BigInt(5));  // 101

  num.clearBit(0);
  EXPECT_EQ(num, BigInt(4));  // 100
}

TEST_F(BigIntBitwiseTest, ClearBitAlreadyClear) {
  BigInt num(5);  // 101
  num.clearBit(1);
  EXPECT_EQ(num, BigInt(5));  // Still 101
}

TEST_F(BigIntBitwiseTest, ClearBitOutOfRange) {
  BigInt num(255);
  num.clearBit(1000);  // Should be no-op
  EXPECT_EQ(num, BigInt(255));
}

TEST_F(BigIntBitwiseTest, ClearBitToZero) {
  BigInt num(1);
  num.clearBit(0);
  EXPECT_EQ(num, BigInt(0));
}

TEST_F(BigIntBitwiseTest, ClearBitChaining) {
  BigInt num(7);  // 111
  num.clearBit(0).clearBit(1).clearBit(2);
  EXPECT_EQ(num, BigInt(0));
}

// ============================================================================
// flipBit Tests
// ============================================================================

TEST_F(BigIntBitwiseTest, FlipBitBasic) {
  BigInt num(5);  // 101
  num.flipBit(1);
  EXPECT_EQ(num, BigInt(7));  // 111

  num.flipBit(0);
  EXPECT_EQ(num, BigInt(6));  // 110
}

TEST_F(BigIntBitwiseTest, FlipBitTwice) {
  BigInt num(42);
  BigInt original = num;
  num.flipBit(3).flipBit(3);
  EXPECT_EQ(num, original);
}

TEST_F(BigIntBitwiseTest, FlipBitExtends) {
  BigInt num(1);
  num.flipBit(100);
  EXPECT_TRUE(num.testBit(0));
  EXPECT_TRUE(num.testBit(100));
}

TEST_F(BigIntBitwiseTest, FlipBitChaining) {
  BigInt num(0);
  num.flipBit(0).flipBit(1).flipBit(2);
  EXPECT_EQ(num, BigInt(7));
}

// ============================================================================
// divmod Tests
// ============================================================================

TEST_F(BigIntBitwiseTest, DivmodBasic) {
  BigInt a(17);
  BigInt b(5);
  auto [q, r] = divmod(a, b);
  EXPECT_EQ(q, BigInt(3));
  EXPECT_EQ(r, BigInt(2));
}

TEST_F(BigIntBitwiseTest, DivmodExact) {
  BigInt a(100);
  BigInt b(10);
  auto [q, r] = divmod(a, b);
  EXPECT_EQ(q, BigInt(10));
  EXPECT_EQ(r, BigInt(0));
}

TEST_F(BigIntBitwiseTest, DivmodLarge) {
  BigInt a("123456789012345678901234567890", 10);
  BigInt b("9876543210", 10);
  auto [q, r] = divmod(a, b);

  // Verify: a == q * b + r
  EXPECT_EQ(a, q * b + r);
  // Verify: 0 <= r < b
  EXPECT_TRUE(r >= BigInt(0));
  EXPECT_TRUE(r < b);
}

TEST_F(BigIntBitwiseTest, DivmodNegative) {
  BigInt a(-17);
  BigInt b(5);
  auto [q, r] = divmod(a, b);
  // -17 / 5 = -3, -17 % 5 = -2
  EXPECT_EQ(q, BigInt(-3));
  EXPECT_EQ(r, BigInt(-2));
}

TEST_F(BigIntBitwiseTest, DivmodByZeroThrows) {
  BigInt a(100);
  BigInt zero(0);
  EXPECT_THROW(divmod(a, zero), std::domain_error);
}

TEST_F(BigIntBitwiseTest, DivmodConsistency) {
  // Test that divmod is consistent with / and %
  BigInt a("999999999999999999999", 10);
  BigInt b("123456789", 10);
  auto [q, r] = divmod(a, b);
  EXPECT_EQ(q, a / b);
  EXPECT_EQ(r, a % b);
}
