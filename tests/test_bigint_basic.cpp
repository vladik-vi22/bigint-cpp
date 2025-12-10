#include <gtest/gtest.h>
#include <bigint/BigInt.hpp>
#include <tuple>

using namespace bigint;

class BigIntBasicTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(BigIntBasicTest, DefaultConstructor) {
    BigInt num;
    EXPECT_EQ(num.toStdString(10), "0");
    EXPECT_TRUE(num.isZero());
}

// ============================================================================
// Parameterized Tests for Base Conversions
// ============================================================================

struct BaseConversionTestCase {
    std::string input;
    int base;
    std::string expectedDecimal;
};

class BigIntBaseConversionTest
    : public ::testing::TestWithParam<BaseConversionTestCase> {};

TEST_P(BigIntBaseConversionTest, StringConstructorWithBase) {
    const auto& tc = GetParam();
    BigInt num(tc.input, tc.base);
    EXPECT_EQ(num.toStdString(10), tc.expectedDecimal);
}

INSTANTIATE_TEST_SUITE_P(
    BaseConversions,
    BigIntBaseConversionTest,
    ::testing::Values(
        // Binary (base 2)
        BaseConversionTestCase{"0", 2, "0"},
        BaseConversionTestCase{"1", 2, "1"},
        BaseConversionTestCase{"10", 2, "2"},
        BaseConversionTestCase{"11111111", 2, "255"},
        BaseConversionTestCase{"100000000", 2, "256"},
        // Note: Octal (base 8) is not supported by this library
        // Decimal (base 10)
        BaseConversionTestCase{"0", 10, "0"},
        BaseConversionTestCase{"255", 10, "255"},
        BaseConversionTestCase{"12345678901234567890", 10, "12345678901234567890"},
        // Hexadecimal (base 16)
        BaseConversionTestCase{"0", 16, "0"},
        BaseConversionTestCase{"F", 16, "15"},
        BaseConversionTestCase{"FF", 16, "255"},
        BaseConversionTestCase{"100", 16, "256"},
        BaseConversionTestCase{"FFFFFFFF", 16, "4294967295"},
        BaseConversionTestCase{"FFFFFFFFFFFFFFFF", 16, "18446744073709551615"}
    )
);

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
    EXPECT_EQ(num.bitLength(), 8);
}

// ============================================================================
// isZero() Tests
// ============================================================================

TEST_F(BigIntBasicTest, IsZeroDefault) {
    BigInt num;
    EXPECT_TRUE(num.isZero());
}

TEST_F(BigIntBasicTest, IsZeroFromString) {
    BigInt zero("0", 10);
    EXPECT_TRUE(zero.isZero());

    BigInt nonZero("1", 10);
    EXPECT_FALSE(nonZero.isZero());
}

TEST_F(BigIntBasicTest, IsZeroAfterOperations) {
    BigInt a("100", 10);
    BigInt b("100", 10);
    BigInt result = a - b;
    EXPECT_TRUE(result.isZero());

    BigInt c("0", 10);
    BigInt d("5", 10);
    BigInt product = c * d;
    EXPECT_TRUE(product.isZero());
}

TEST_F(BigIntBasicTest, IsZeroNegativeZero) {
    // Ensure -0 is treated as zero
    BigInt negZero("-0", 10);
    EXPECT_TRUE(negZero.isZero());
}

// ============================================================================
// Edge Case Tests
// ============================================================================

TEST_F(BigIntBasicTest, ZeroOperations) {
    BigInt zero;
    BigInt num("12345", 10);

    // Addition with zero
    EXPECT_EQ((num + zero).toStdString(10), "12345");
    EXPECT_EQ((zero + num).toStdString(10), "12345");

    // Subtraction with zero
    EXPECT_EQ((num - zero).toStdString(10), "12345");

    // Multiplication with zero
    EXPECT_TRUE((num * zero).isZero());
    EXPECT_TRUE((zero * num).isZero());
}

TEST_F(BigIntBasicTest, OneOperations) {
    BigInt one = constants::One;
    BigInt num("12345", 10);

    // Multiplication by one
    EXPECT_EQ((num * one).toStdString(10), "12345");
    EXPECT_EQ((one * num).toStdString(10), "12345");

    // Division by one
    EXPECT_EQ((num / one).toStdString(10), "12345");
}

TEST_F(BigIntBasicTest, SignHandling) {
    BigInt pos("100", 10);
    BigInt neg("-100", 10);

    // Positive * Negative = Negative
    EXPECT_TRUE((pos * neg).isNegative());

    // Negative * Negative = Positive
    EXPECT_TRUE((neg * neg).isPositive());

    // Negative + Negative = Negative
    EXPECT_TRUE((neg + neg).isNegative());
}

TEST_F(BigIntBasicTest, BoundaryValues) {
    // Max uint32_t
    BigInt maxU32("4294967295", 10);
    EXPECT_EQ(maxU32.toStdString(16), "ffffffff");  // Library uses lowercase hex

    // Max uint64_t
    BigInt maxU64("18446744073709551615", 10);
    EXPECT_EQ(maxU64.toStdString(16), "ffffffffffffffff");

    // Just over uint64_t
    BigInt overU64("18446744073709551616", 10);
    EXPECT_EQ(overU64.toStdString(16), "10000000000000000");
}

TEST_F(BigIntBasicTest, BitLengthEdgeCases) {
    // Note: Library returns 1 for zero's bit length (represents single 0 bit)
    EXPECT_EQ(BigInt("0", 10).bitLength(), 1);
    EXPECT_EQ(BigInt("1", 10).bitLength(), 1);
    EXPECT_EQ(BigInt("2", 10).bitLength(), 2);
    EXPECT_EQ(BigInt("3", 10).bitLength(), 2);
    EXPECT_EQ(BigInt("4", 10).bitLength(), 3);
    EXPECT_EQ(BigInt("128", 10).bitLength(), 8);
    EXPECT_EQ(BigInt("256", 10).bitLength(), 9);
}

TEST_F(BigIntBasicTest, ByteLengthEdgeCases) {
    // Note: Library returns 1 for zero's byte length
    EXPECT_EQ(BigInt("0", 10).byteLength(), 1);
    EXPECT_EQ(BigInt("1", 10).byteLength(), 1);
    EXPECT_EQ(BigInt("255", 10).byteLength(), 1);
    EXPECT_EQ(BigInt("256", 10).byteLength(), 2);
    EXPECT_EQ(BigInt("65535", 10).byteLength(), 2);
    EXPECT_EQ(BigInt("65536", 10).byteLength(), 3);
}

