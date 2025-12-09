#include <gtest/gtest.h>
#include <bigint/BigInt.hpp>

using namespace bigint;

class BigIntBasicTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(BigIntBasicTest, DefaultConstructor) {
    BigInt num;
    EXPECT_EQ(num.toStdString(10), "0");
}

TEST_F(BigIntBasicTest, StringConstructorDecimal) {
    BigInt num("12345678901234567890", 10);
    EXPECT_EQ(num.toStdString(10), "12345678901234567890");
}

TEST_F(BigIntBasicTest, StringConstructorHex) {
    BigInt num("FF", 16);
    EXPECT_EQ(num.toStdString(10), "255");
}

TEST_F(BigIntBasicTest, StringConstructorBinary) {
    BigInt num("11111111", 2);
    EXPECT_EQ(num.toStdString(10), "255");
}

TEST_F(BigIntBasicTest, CopyConstructor) {
    BigInt original("999999999999999999", 10);
    BigInt copy(original);
    EXPECT_EQ(copy.toStdString(10), original.toStdString(10));
}

TEST_F(BigIntBasicTest, Uint64Constructor) {
    BigInt num(static_cast<uint64_t>(18446744073709551615ULL));
    EXPECT_EQ(num.toStdString(10), "18446744073709551615");
}

TEST_F(BigIntBasicTest, NegativeNumber) {
    BigInt num("-12345", 10);
    EXPECT_TRUE(num.isNegative());
    EXPECT_FALSE(num.isPositive());
}

TEST_F(BigIntBasicTest, IsEvenOdd) {
    BigInt even("100", 10);
    BigInt odd("101", 10);
    
    EXPECT_TRUE(even.isEven());
    EXPECT_FALSE(even.isOdd());
    EXPECT_TRUE(odd.isOdd());
    EXPECT_FALSE(odd.isEven());
}

TEST_F(BigIntBasicTest, BitLength) {
    BigInt num("255", 10);  // 11111111 in binary = 8 bits
    EXPECT_EQ(num.bitLenght(), 8);
}

