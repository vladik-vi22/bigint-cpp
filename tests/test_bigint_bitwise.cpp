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
  EXPECT_TRUE((zero << static_cast<size_t>(10)).isZero());
  EXPECT_TRUE((zero >> static_cast<size_t>(10)).isZero());
}

TEST_F(BigIntBitwiseTest, RightShiftToZero) {
  BigInt num("255", 10);
  BigInt result = num >> static_cast<size_t>(100);  // Shift more than bit length
  EXPECT_TRUE(result.isZero());
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

  EXPECT_TRUE((num & zero).isZero());
  EXPECT_EQ((num | zero).toStdString(10), "12345");
  EXPECT_EQ((num ^ zero).toStdString(10), "12345");
}

TEST_F(BigIntBitwiseTest, BitwiseWithSelf) {
  BigInt num("12345", 10);

  EXPECT_EQ((num & num).toStdString(10), "12345");
  EXPECT_EQ((num | num).toStdString(10), "12345");
  EXPECT_TRUE((num ^ num).isZero());
}
