#pragma once

#include <array>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace bigint {

/**
 * @brief Arbitrary precision signed integer.
 *
 * Internal representation uses a little-endian vector of 32-bit words
 * with sign-magnitude encoding. Supports arithmetic, bitwise, and
 * comparison operations with automatic memory management.
 */
class BigInt
{
public:
    /// @name Base Constants
    /// @{
    static constexpr uint8_t baseBinary = 2;       ///< Binary base for string I/O
    static constexpr uint8_t baseDecimal = 10;     ///< Decimal base for string I/O
    static constexpr uint8_t baseHexadecimal = 16; ///< Hexadecimal base for string I/O
    static constexpr uint64_t basisCalcSys = 1ULL << 32;  ///< Internal calculation base (2^32)
    /// @}

    static inline uint8_t baseInput = baseDecimal;  ///< Default base for input parsing
    static inline uint8_t baseOutput = baseDecimal; ///< Default base for output formatting

private:
    bool positive;                      ///< Sign flag (true = positive or zero)
    std::vector<uint32_t> vectorUint32_t; ///< Little-endian 32-bit word storage

public:
    /// @name Constructors
    /// @{

    /// @brief Default constructor. Initializes to zero.
    BigInt();

    /// @brief Copy constructor.
    BigInt(const BigInt& bigInt);

    /// @brief Construct from string representation.
    /// @param bigIntString Number as string (may include leading '-' for negative)
    /// @param base Numeric base (2-16, default: baseInput)
    explicit BigInt(std::string bigIntString, uint8_t base = baseInput);

    /// @brief Construct from fixed-size array of 32-bit words (big-endian order).
    template<size_t size>
    explicit BigInt(const std::array<uint32_t, size>& bigIntArrayUint32_t, bool isPositive = true)
    {
        vectorUint32_t = std::vector<uint32_t>(bigIntArrayUint32_t.crbegin(), bigIntArrayUint32_t.crend());
        positive = isPositive;
    }

    /// @brief Construct from vector of 32-bit words (big-endian order).
    explicit BigInt(const std::vector<uint32_t>& bigIntVectorUint32_t, bool isPositive = true);

    /// @brief Construct from vector of 16-bit words (big-endian order).
    explicit BigInt(const std::vector<uint16_t>& bigIntVectorUint16_t, bool isPositive = true);

    /// @brief Construct from vector of bytes (big-endian order).
    explicit BigInt(const std::vector<uint8_t>& bigIntVectorUint8_t, bool isPositive = true);

    /// @brief Construct from vector of bits (big-endian order).
    explicit BigInt(const std::vector<bool>& bigIntVectorBool, bool isPositive = true);

    /// @brief Construct from 64-bit unsigned integer.
    explicit BigInt(uint64_t bigIntUint64_t, bool isPositive = true);

    /// @brief Construct from 32-bit unsigned integer.
    explicit BigInt(uint32_t bigIntUint64_t, bool isPositive = true);

    /// @brief Construct from 64-bit signed integer.
    explicit BigInt(int64_t bigIntInt64_t);

    /// @brief Construct from 32-bit signed integer.
    explicit BigInt(int32_t bigIntInt32_t);

    ~BigInt() = default;
    /// @}

    /// @name Assignment
    /// @{
    BigInt& operator = (const BigInt& equal);
    /// @}

    /// @name Arithmetic Operators
    /// @{
    BigInt operator +() const;                        ///< Unary plus
    BigInt operator + (const BigInt& addend) const;   ///< Addition
    BigInt& operator += (const BigInt& addend);       ///< Addition assignment
    BigInt& operator ++();                            ///< Pre-increment
    BigInt operator ++(int);                          ///< Post-increment

    BigInt operator -() const;                            ///< Unary minus (negation)
    BigInt operator - (const BigInt& subtrahend) const;   ///< Subtraction
    BigInt& operator -= (const BigInt& subtrahend);       ///< Subtraction assignment
    BigInt& operator --();                                ///< Pre-decrement
    BigInt operator --(int);                              ///< Post-decrement

    BigInt operator * (const BigInt& multiplier) const;   ///< Multiplication
    BigInt& operator *= (const BigInt& multiplier);       ///< Multiplication assignment

    BigInt operator / (const BigInt& divisor) const;      ///< Integer division
    BigInt& operator /= (const BigInt& divisor);          ///< Division assignment

    BigInt operator % (const BigInt& divisor) const;      ///< Modulo (remainder)
    BigInt& operator %= (const BigInt& divisor);          ///< Modulo assignment
    /// @}

    /// @name Mathematical Functions
    /// @{

    /// @brief Compute base raised to exponent.
    friend BigInt pow(const BigInt& base, const BigInt& exponent);

    /// @brief Compute floor of log base 2.
    friend size_t log2(const BigInt& antilogarithm);

    /// @brief Compute (base^exponent) mod divisor efficiently.
    friend BigInt powmod(const BigInt& base, const BigInt& exponent, const BigInt& divisor);

    /// @brief Compute modular multiplicative inverse.
    friend BigInt inversemod(BigInt dividend, const BigInt& divisor);

    /// @brief Check if dividend1 ≡ dividend2 (mod divisor).
    friend bool congruencemod(const BigInt& dividend1, const BigInt& dividend2, BigInt divisor);

    /// @brief Check if two numbers are coprime (gcd == 1).
    friend bool isCoprime(const BigInt& bigInt1, const BigInt& bigInt2);

    /// @brief Compute Jacobi symbol (a/n).
    friend int8_t symbolJacobi(BigInt bigInt1, BigInt bigInt2);
    /// @}

    /// @name Bitwise Operators
    /// @{
    BigInt operator ~() const;                                  ///< Bitwise NOT
    BigInt operator & (const BigInt& rightBitwiseAND) const;    ///< Bitwise AND
    BigInt& operator &= (const BigInt& rightBitwiseAND);        ///< Bitwise AND assignment
    BigInt operator | (const BigInt& rightBitwiseOR) const;     ///< Bitwise OR
    BigInt& operator |= (const BigInt& rightBitwiseOR);         ///< Bitwise OR assignment
    BigInt operator ^ (const BigInt& rightBitwiseXOR) const;    ///< Bitwise XOR
    BigInt& operator ^= (const BigInt& rightBitwiseXOR);        ///< Bitwise XOR assignment

    BigInt operator << (size_t shift) const;    ///< Left shift
    BigInt& operator <<= (size_t shift);        ///< Left shift assignment
    BigInt operator >> (size_t shift) const;    ///< Right shift
    BigInt& operator >>= (size_t shift);        ///< Right shift assignment

    /// @brief Circular left shift within current bit length.
    [[nodiscard]] BigInt leftCircularShift(size_t shift) const;

    /// @brief Circular right shift within current bit length.
    [[nodiscard]] BigInt rightCircularShift(size_t shift) const;
    /// @}

    /// @name Logical Operators
    /// @{
    bool operator !() const;                        ///< Logical NOT (true if zero)
    bool operator && (const BigInt& rightAND) const;  ///< Logical AND
    bool operator || (const BigInt& rightOR) const;   ///< Logical OR
    /// @}

    /// @name Comparison Operators
    /// @{
    bool operator == (const BigInt& rightComparable) const;
    bool operator != (const BigInt& rightComparable) const;
    bool operator < (const BigInt& rightComparable) const;
    bool operator > (const BigInt& rightComparable) const;
    bool operator <= (const BigInt& rightComparable) const;
    bool operator >= (const BigInt& rightComparable) const;
    /// @}

    /// @name Utility Functions
    /// @{
    friend BigInt abs(const BigInt& bigInt);                      ///< Absolute value
    friend BigInt gcd(BigInt bigInt1, BigInt bigInt2);            ///< Greatest common divisor
    friend BigInt lcm(const BigInt& bigInt1, const BigInt& bigInt2); ///< Least common multiple
    friend const BigInt& max(const BigInt& bigInt1, const BigInt& bigInt2);
    friend const BigInt& min(const BigInt& bigInt1, const BigInt& bigInt2);
    /// @}

    /// @name Stream I/O
    /// @{
    friend std::ostream& operator << (std::ostream& out, const BigInt& bigInt);
    friend std::istream& operator >> (std::istream& in, BigInt& bigInt);
    /// @}

    /// @name Conversion
    /// @{

    /// @brief Convert to string in specified base (2-16).
    [[nodiscard]] std::string toStdString(int base = baseOutput) const;

    /// @brief Convert to vector of 32-bit words (big-endian order).
    [[nodiscard]] std::vector<uint32_t> toStdVectorUint32_t() const;

    /// @brief Convert to vector of bytes (big-endian order).
    [[nodiscard]] std::vector<uint8_t> toStdVectorUint8_t() const;

    [[nodiscard]] explicit operator uint64_t() const;  ///< Convert to uint64_t (truncates)
    [[nodiscard]] explicit operator uint32_t() const;  ///< Convert to uint32_t (truncates)
    [[nodiscard]] explicit operator uint16_t() const;  ///< Convert to uint16_t (truncates)
    [[nodiscard]] explicit operator uint8_t() const;   ///< Convert to uint8_t (truncates)
    [[nodiscard]] explicit operator bool() const;      ///< True if non-zero
    /// @}

    /// @name Queries
    /// @{
    [[nodiscard]] size_t bitLength() const;   ///< Number of bits (minimum 1 for zero)
    [[nodiscard]] size_t byteLength() const;  ///< Number of bytes needed to represent value
    [[nodiscard]] bool isZero() const;        ///< True if value is zero
    [[nodiscard]] bool isEven() const;        ///< True if value is even
    [[nodiscard]] bool isOdd() const;         ///< True if value is odd
    [[nodiscard]] bool isPositive() const;    ///< True if value >= 0
    [[nodiscard]] bool isNegative() const;    ///< True if value < 0
    /// @}

private:
    BigInt operator * (uint32_t multiplier) const;
    BigInt& operator *= (uint32_t multiplier);
    std::pair<BigInt, BigInt> DivMod(const BigInt& divisor) const;
    friend BigInt BarrettReduction(const BigInt& dividend, const BigInt& divisor, const BigInt& mu);
    void alignTo(BigInt& aligned);
    void deleteZeroHighOrderDigit();
    BigInt shiftDigitsToHigh(size_t shift) const;
    BigInt shiftDigitsToLow(size_t shift) const;
    BigInt toBigIntDec() const;
};

/// @brief Common BigInt constants for convenience.
namespace constants {
    inline const BigInt Zero{0};   ///< Constant 0
    inline const BigInt One{1};    ///< Constant 1
    inline const BigInt Two{2};    ///< Constant 2
    inline const BigInt Three{3};  ///< Constant 3
    inline const BigInt Four{4};   ///< Constant 4
    inline const BigInt Five{5};   ///< Constant 5
    inline const BigInt Eight{8};  ///< Constant 8
}

} // namespace bigint

