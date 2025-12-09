#include <gtest/gtest.h>
#include <bigint/BigInt.hpp>

using namespace bigint;

class BigIntBitwiseTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(BigIntBitwiseTest, BitwiseAND) {
    BigInt a("15", 10);   // 1111
    BigInt b("9", 10);    // 1001
    BigInt result = a & b;
    EXPECT_EQ(result.toStdString(10), "9");  // 1001
}

TEST_F(BigIntBitwiseTest, BitwiseOR) {
    BigInt a("12", 10);   // 1100
    BigInt b("10", 10);   // 1010
    BigInt result = a | b;
    EXPECT_EQ(result.toStdString(10), "14");  // 1110
}

TEST_F(BigIntBitwiseTest, BitwiseXOR) {
    BigInt a("12", 10);   // 1100
    BigInt b("10", 10);   // 1010
    BigInt result = a ^ b;
    EXPECT_EQ(result.toStdString(10), "6");   // 0110
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

TEST_F(BigIntBitwiseTest, Comparison) {
    BigInt a("100", 10);
    BigInt b("200", 10);
    BigInt c("100", 10);
    
    EXPECT_TRUE(a < b);
    EXPECT_TRUE(b > a);
    EXPECT_TRUE(a == c);
    EXPECT_TRUE(a != b);
    EXPECT_TRUE(a <= c);
    EXPECT_TRUE(a >= c);
    EXPECT_TRUE(a <= b);
}

