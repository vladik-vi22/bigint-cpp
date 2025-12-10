#pragma once

#include <array>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace bigint {

/// Arbitrary precision signed integer.
/// Internal representation: little-endian vector of 32-bit words, sign-magnitude.
class BigInt
{
public:
    // --- Base constants ---
    static constexpr uint8_t baseBinary = 2;
    static constexpr uint8_t baseDecimal = 10;
    static constexpr uint8_t baseHexadecimal = 16;
    static constexpr uint64_t basisCalcSys = 1ULL << 32;  // 2^32

    static inline uint8_t baseInput = baseDecimal;
    static inline uint8_t baseOutput = baseDecimal;

private:
    bool positive;
    std::vector<uint32_t> vectorUint32_t;

public:
    BigInt();
    BigInt(const BigInt& bigInt);
    explicit BigInt(std::string bigIntString, uint8_t base = baseInput);
    template<size_t size>
    explicit BigInt(const std::array<uint32_t, size>& bigIntArrayUint32_t, bool isPositive = true)
    {
        vectorUint32_t = std::vector<uint32_t>(bigIntArrayUint32_t.crbegin(), bigIntArrayUint32_t.crend());
        positive = isPositive;
    }
    explicit BigInt(const std::vector<uint32_t>& bigIntVectorUint32_t, bool isPositive = true);
    explicit BigInt(const std::vector<uint16_t>& bigIntVectorUint16_t, bool isPositive = true);
    explicit BigInt(const std::vector<uint8_t>& bigIntVectorUint8_t, bool isPositive = true);
    explicit BigInt(const std::vector<bool>& bigIntVectorBool, bool isPositive = true);
    explicit BigInt(uint64_t bigIntUint64_t, bool isPositive = true);
    explicit BigInt(uint32_t bigIntUint64_t, bool isPositive = true);
    explicit BigInt(int64_t bigIntInt64_t);
    explicit BigInt(int32_t bigIntInt32_t);
    ~BigInt() = default;

    BigInt& operator = (const BigInt& equal); // simple assignment

    BigInt operator +() const; // unary plus
    BigInt operator + (const BigInt& addend) const; // addition
    BigInt& operator += (const BigInt& addend); // addition assignment
    BigInt& operator ++(); // pre-increment
    BigInt operator ++(int); // post-increment

    BigInt operator -() const; // unary minus
    BigInt operator - (const BigInt& subtrahend) const; // subtraction
    BigInt& operator -= (const BigInt& subtrahend); // subtraction assignment
    BigInt& operator --(); // pre-decrement
    BigInt operator --(int); // post-decrement

    BigInt operator * (const BigInt& multiplier) const; // multiplication
    BigInt& operator *= (const BigInt& multiplier); // multiplication assignment

    BigInt operator / (const BigInt& divisor) const; // division
    BigInt& operator /= (const BigInt& divisor); // division assignment

    BigInt operator % (const BigInt& divisor) const; // modulo
    BigInt& operator %= (const BigInt& divisor); // modulo assignment

    friend BigInt pow(const BigInt& base, const BigInt& exponent); // power

    friend size_t log2(const BigInt& antilogarithm); // logarithm to the base 2

    friend BigInt powmod(const BigInt& base, const BigInt& exponent, const BigInt& divisor);
    friend BigInt inversemod(BigInt dividend, const BigInt& divisor);
    friend bool congruencemod(const BigInt& dividend1, const BigInt& dividend2, BigInt divisor);
    friend bool isCoprime(const BigInt& bigInt1, const BigInt& bigInt2);

    friend int8_t symbolJacobi(BigInt bigInt1, BigInt bigInt2);

    BigInt operator ~() const; // bitwise NOT
    BigInt operator & (const BigInt& rightBitwiseAND) const; // bitwise AND
    BigInt& operator &= (const BigInt& rightBitwiseAND); // bitwise AND assignment
    BigInt operator | (const BigInt& rightBitwiseOR) const; // bitwise OR
    BigInt& operator |= (const BigInt& rightBitwiseOR); // bitwise OR assignment
    BigInt operator ^ (const BigInt& rightBitwiseXOR) const; // bitwise XOR
    BigInt& operator ^= (const BigInt& rightBitwiseXOR); // bitwise XOR assignment

    BigInt operator << (size_t shift) const; // bitwise left shift
    BigInt& operator <<= (size_t shift); // bitwise left shift assignment
    BigInt operator >> (size_t shift) const; // bitwise right shift
    BigInt& operator >>= (size_t shift); // bitwise right shift assignment

    [[nodiscard]] BigInt leftCircularShift(size_t shift) const; // bitwise left circular shift
    [[nodiscard]] BigInt rightCircularShift(size_t shift) const; // bitwise right circular shift;

    bool operator !() const; // negation
    bool operator && (const BigInt& rightAND) const; // AND
    bool operator || (const BigInt& rightOR) const; // inclusive OR

    bool operator == (const BigInt& rightComparable) const; // equal to
    bool operator != (const BigInt& rightComparable) const; // not equal to
    bool operator < (const BigInt& rightComparable) const; // less than
    bool operator > (const BigInt& rightComparable) const; // greater than
    bool operator <= (const BigInt& rightComparable) const; // less than or equal to
    bool operator >= (const BigInt& rightComparable) const; // greater than or equal to

    friend BigInt abs(const BigInt& bigInt); // absolute value
    friend BigInt gcd(BigInt bigInt1, BigInt bigInt2); // greatest common divisor
    friend BigInt lcm(const BigInt& bigInt1, const BigInt& bigInt2); // least common multiple

    friend const BigInt& max(const BigInt& bigInt1, const BigInt& bigInt2);
    friend const BigInt& min(const BigInt& bigInt1, const BigInt& bigInt2);

    friend std::ostream& operator << (std::ostream& out, const BigInt& bigInt);
    friend std::istream& operator >> (std::istream& in, BigInt& bigInt);

    // --- Conversion ---
    [[nodiscard]] std::string toStdString(int base = baseOutput) const;
    [[nodiscard]] std::vector<uint32_t> toStdVectorUint32_t() const;
    [[nodiscard]] std::vector<uint8_t> toStdVectorUint8_t() const;
    [[nodiscard]] explicit operator uint64_t() const;
    [[nodiscard]] explicit operator uint32_t() const;
    [[nodiscard]] explicit operator uint16_t() const;
    [[nodiscard]] explicit operator uint8_t() const;
    [[nodiscard]] explicit operator bool() const;

    // --- Queries ---
    [[nodiscard]] size_t bitLength() const;
    [[nodiscard]] size_t byteLength() const;
    [[nodiscard]] bool isZero() const;
    [[nodiscard]] bool isEven() const;
    [[nodiscard]] bool isOdd() const;
    [[nodiscard]] bool isPositive() const;
    [[nodiscard]] bool isNegative() const;

private:
    BigInt operator * (uint32_t multiplier) const; // multiplication
    BigInt& operator *= (uint32_t multiplier); // multiplication assignment
    std::pair<BigInt, BigInt> DivMod(const BigInt& divisor) const;
    friend BigInt BarrettReduction(const BigInt& dividend, const BigInt& divisor, const BigInt& mu);
    void alignTo(BigInt& aligned);
    void deleteZeroHighOrderDigit();
    BigInt shiftDigitsToHigh(size_t shift) const;
    BigInt shiftDigitsToLow(size_t shift) const;
    BigInt toBigIntDec() const;
};

/// Common constants for convenience.
namespace constants {
    inline const BigInt Zero{0};
    inline const BigInt One{1};
    inline const BigInt Two{2};
    inline const BigInt Three{3};
    inline const BigInt Four{4};
    inline const BigInt Five{5};
    inline const BigInt Eight{8};
}

std::string strDec2strBin(std::string strDec);

} // namespace bigint

