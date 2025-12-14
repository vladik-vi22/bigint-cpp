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
