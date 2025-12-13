/**
 * @file BigInt.hpp
 * @brief Arbitrary precision integer class.
 */
#pragma once

#include <bigint/Constants.hpp>

#include <array>
#include <compare>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

namespace bigint {

/**
 * @brief Arbitrary precision signed integer.
 *
 * Internal representation uses a little-endian vector of 32-bit words
 * with sign-magnitude encoding. Supports arithmetic, bitwise, and
 * comparison operations with automatic memory management.
 */
class BigInt {
 public:
  /// @name Constructors
  /// @{

  /// @brief Default constructor. Initializes to zero.
  BigInt();

  /// @brief Copy constructor.
  BigInt(const BigInt& other);

  /// @brief Move constructor.
  BigInt(BigInt&& other) noexcept;

  /// @brief Construct from string representation.
  /// @param str Number as string (may include leading '-' for negative)
  /// @param base Numeric base (2, 10, or 16; default: kDefaultInputBase)
  explicit BigInt(std::string str, uint8_t base = kDefaultInputBase);

  /// @brief Construct from fixed-size array of 32-bit words (big-endian order).
  template <size_t N>
  explicit BigInt(const std::array<uint32_t, N>& arr, bool is_positive = true) {
    digits_ = std::vector<uint32_t>(arr.crbegin(), arr.crend());
    positive_ = is_positive;
  }

  /// @brief Construct from vector of 32-bit words (big-endian order).
  explicit BigInt(const std::vector<uint32_t>& vec, bool is_positive = true);

  /// @brief Construct from vector of 16-bit words (big-endian order).
  explicit BigInt(const std::vector<uint16_t>& vec, bool is_positive = true);

  /// @brief Construct from vector of bytes (big-endian order).
  explicit BigInt(const std::vector<uint8_t>& vec, bool is_positive = true);

  /// @brief Construct from vector of bits (big-endian order).
  explicit BigInt(const std::vector<bool>& vec, bool is_positive = true);

  /// @brief Construct from 64-bit unsigned integer.
  explicit BigInt(uint64_t value, bool is_positive = true);

  /// @brief Construct from 32-bit unsigned integer.
  explicit BigInt(uint32_t value, bool is_positive = true);

  /// @brief Construct from 64-bit signed integer.
  explicit BigInt(int64_t value);

  /// @brief Construct from 32-bit signed integer.
  explicit BigInt(int32_t value);

  ~BigInt() = default;
  /// @}

  /// @name Assignment
  /// @{
  BigInt& operator=(const BigInt& other);
  BigInt& operator=(BigInt&& other) noexcept;
  /// @}

  /// @name Arithmetic Operators
  /// @{
  [[nodiscard]] BigInt operator+() const;                       ///< Unary plus
  [[nodiscard]] BigInt operator+(const BigInt& addend) const;   ///< Addition
  BigInt& operator+=(const BigInt& addend);       ///< Addition assignment
  BigInt& operator++();                           ///< Pre-increment
  BigInt operator++(int);                         ///< Post-increment

  [[nodiscard]] BigInt operator-() const;                           ///< Unary minus (negation)
  [[nodiscard]] BigInt operator-(const BigInt& subtrahend) const;   ///< Subtraction
  BigInt& operator-=(const BigInt& subtrahend);       ///< Subtraction assignment
  BigInt& operator--();                               ///< Pre-decrement
  BigInt operator--(int);                             ///< Post-decrement

  [[nodiscard]] BigInt operator*(const BigInt& multiplier) const;   ///< Multiplication
  BigInt& operator*=(const BigInt& multiplier);       ///< Multiplication assignment

  [[nodiscard]] BigInt operator/(const BigInt& divisor) const;      ///< Integer division
  BigInt& operator/=(const BigInt& divisor);          ///< Division assignment

  [[nodiscard]] BigInt operator%(const BigInt& divisor) const;      ///< Modulo (remainder)
  BigInt& operator%=(const BigInt& divisor);          ///< Modulo assignment
  /// @}

  /// @name Mathematical Functions
  /// @{

  /// @brief Compute base raised to exponent.
  friend BigInt pow(const BigInt& base, const BigInt& exponent);

  /// @brief Compute floor of log base 2.
  friend size_t log2(const BigInt& antilogarithm);

  /// @brief Compute (base^exponent) mod divisor efficiently.
  friend BigInt powmod(const BigInt& base, const BigInt& exponent,
                       const BigInt& divisor);

  /// @brief Compute modular multiplicative inverse.
  friend BigInt inversemod(BigInt dividend, const BigInt& divisor);

  /// @brief Check if dividend1 = dividend2 (mod divisor).
  friend bool congruencemod(const BigInt& dividend1, const BigInt& dividend2,
                            BigInt divisor);

  /// @brief Check if two numbers are coprime (gcd == 1).
  friend bool isCoprime(const BigInt& a, const BigInt& b);

  /// @brief Compute Jacobi symbol (a/n).
  friend int8_t symbolJacobi(BigInt a, BigInt b);
  /// @}

  /// @name Bitwise Operators
  /// @{
  [[nodiscard]] BigInt operator~() const;                   ///< Bitwise NOT
  [[nodiscard]] BigInt operator&(const BigInt& rhs) const;  ///< Bitwise AND
  BigInt& operator&=(const BigInt& rhs);      ///< Bitwise AND assignment
  [[nodiscard]] BigInt operator|(const BigInt& rhs) const;  ///< Bitwise OR
  BigInt& operator|=(const BigInt& rhs);      ///< Bitwise OR assignment
  [[nodiscard]] BigInt operator^(const BigInt& rhs) const;  ///< Bitwise XOR
  BigInt& operator^=(const BigInt& rhs);      ///< Bitwise XOR assignment

  [[nodiscard]] BigInt operator<<(size_t shift) const;  ///< Left shift
  BigInt& operator<<=(size_t shift);      ///< Left shift assignment
  [[nodiscard]] BigInt operator>>(size_t shift) const;  ///< Right shift
  BigInt& operator>>=(size_t shift);      ///< Right shift assignment

  /// @brief Circular left shift within current bit length.
  [[nodiscard]] BigInt leftCircularShift(size_t shift) const;

  /// @brief Circular right shift within current bit length.
  [[nodiscard]] BigInt rightCircularShift(size_t shift) const;
  /// @}

  /// @name Logical Operators
  /// @{
  bool operator!() const noexcept;                    ///< Logical NOT (true if zero)
  bool operator&&(const BigInt& rhs) const noexcept;  ///< Logical AND
  bool operator||(const BigInt& rhs) const noexcept;  ///< Logical OR
  /// @}

  /// @name Comparison Operators
  /// @{
  [[nodiscard]] std::strong_ordering operator<=>(const BigInt& rhs) const noexcept;
  [[nodiscard]] bool operator==(const BigInt& rhs) const noexcept;
  /// @}

  /// @name Utility Functions
  /// @{
  friend BigInt abs(const BigInt& value);               ///< Absolute value
  friend BigInt gcd(BigInt a, BigInt b);                ///< Greatest common divisor
  friend BigInt lcm(const BigInt& a, const BigInt& b);  ///< Least common multiple
  friend const BigInt& max(const BigInt& a, const BigInt& b) noexcept;  ///< Maximum
  friend const BigInt& min(const BigInt& a, const BigInt& b) noexcept;  ///< Minimum
  friend void swap(BigInt& lhs, BigInt& rhs) noexcept;  ///< Swap two values
  /// @}

  /// @name Stream I/O
  /// @{
  friend std::ostream& operator<<(std::ostream& out, const BigInt& value);
  friend std::istream& operator>>(std::istream& in, BigInt& value);
  /// @}

  /// @name Conversion
  /// @{

  /// @brief Convert to string in specified base (2, 10, or 16).
  [[nodiscard]] std::string toStdString(uint8_t base = kDefaultOutputBase) const;

  /// @brief Convert to vector of 32-bit words (big-endian order).
  [[nodiscard]] std::vector<uint32_t> toStdVectorUint32_t() const;

  /// @brief Convert to vector of bytes (big-endian order).
  [[nodiscard]] std::vector<uint8_t> toStdVectorUint8_t() const;

  [[nodiscard]] explicit operator uint64_t() const noexcept;  ///< Convert to uint64_t (truncates)
  [[nodiscard]] explicit operator uint32_t() const noexcept;  ///< Convert to uint32_t (truncates)
  [[nodiscard]] explicit operator uint16_t() const noexcept;  ///< Convert to uint16_t (truncates)
  [[nodiscard]] explicit operator uint8_t() const noexcept;   ///< Convert to uint8_t (truncates)
  [[nodiscard]] explicit operator bool() const noexcept;      ///< True if non-zero
  /// @}

  /// @name Queries
  /// @{
  [[nodiscard]] size_t bitLength() const noexcept;   ///< Number of bits (minimum 1 for zero)
  [[nodiscard]] size_t byteLength() const noexcept;  ///< Number of bytes needed
  [[nodiscard]] size_t digitCount() const noexcept;  ///< Number of 32-bit digits
  [[nodiscard]] bool isZero() const noexcept;        ///< True if value is zero
  [[nodiscard]] bool isEven() const noexcept;        ///< True if value is even
  [[nodiscard]] bool isOdd() const noexcept;         ///< True if value is odd
  [[nodiscard]] bool isPositive() const noexcept;    ///< True if value >= 0
  [[nodiscard]] bool isNegative() const noexcept;    ///< True if value < 0
  /// @}

 private:
  bool positive_;                   ///< Sign flag (true = positive or zero)
  std::vector<uint32_t> digits_;    ///< Little-endian 32-bit word storage

  /// @brief Compare magnitudes (absolute values) without allocation.
  /// @return -1 if |*this| < |other|, 0 if equal, 1 if |*this| > |other|
  [[nodiscard]] int compareMagnitude(const BigInt& other) const noexcept;

  BigInt operator*(uint32_t multiplier) const;
  BigInt& operator*=(uint32_t multiplier);
  std::pair<BigInt, BigInt> DivMod(const BigInt& divisor) const;
  friend BigInt BarrettReduction(const BigInt& dividend, const BigInt& divisor,
                                 const BigInt& mu);
  void alignTo(BigInt& aligned);
  void deleteZeroHighOrderDigit();
  BigInt shiftDigitsToHigh(size_t shift) const;
  BigInt shiftDigitsToLow(size_t shift) const;
  BigInt toBigIntDec() const;
};

}  // namespace bigint

namespace std {

/// @brief Hash specialization for BigInt to enable use in unordered containers.
template <>
struct hash<bigint::BigInt> {
  size_t operator()(const bigint::BigInt& value) const noexcept {
    // FNV-1a inspired hash combining sign and all digits
    size_t h = 14695981039346656037ULL;  // FNV offset basis
    h ^= static_cast<size_t>(value.isNegative());
    h *= 1099511628211ULL;  // FNV prime
    auto digits = value.toStdVectorUint32_t();
    for (uint32_t digit : digits) {
      h ^= static_cast<size_t>(digit);
      h *= 1099511628211ULL;
    }
    return h;
  }
};

}  // namespace std
