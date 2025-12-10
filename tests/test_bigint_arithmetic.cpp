#include <gtest/gtest.h>
#include <bigint/BigInt.hpp>

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
    BigInt base("123456789012345678901234567890123456789012345678901234567890", 10);
    BigInt exp("65537", 10);  // Common RSA public exponent
    BigInt mod("1000000000000000000000000000000000000000000000000000000000057", 10);
    BigInt result = powmod(base, exp, mod);
    // Result should be less than modulus
    EXPECT_TRUE(result < mod);
    EXPECT_TRUE(result >= constants::Zero);
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
    EXPECT_TRUE(abs(negative).isPositive());
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
    EXPECT_TRUE(isCoprime(large, constants::One));
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
    BigInt inv = inversemod(a, mod);
    EXPECT_EQ(inv.toStdString(10), "0");  // Returns 0 when no inverse
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

