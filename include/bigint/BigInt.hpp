/**
 * @file BigInt.hpp
 * @brief Arbitrary precision integer class.
 */
#pragma once

#include <bigint/Constants.hpp>

#include <array>
#include <compare>
#include <concepts>
#include <cstdint>
#include <iosfwd>
#include <span>
#include <string>
#include <type_traits>
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
  /// @param str Number as string (may include leading '-' for negative).
  ///            Empty string is treated as zero.
  /// @param base Numeric base (2, 10, or 16; default: kDefaultInputBase).
  /// @throws std::invalid_argument if str contains invalid characters for the base.
  explicit BigInt(std::string str, uint8_t base = kDefaultInputBase);

  /// @brief Construct from span of 32-bit words (big-endian order).
  /// @details Accepts any contiguous range: std::vector, std::array, C arrays, etc.
  explicit BigInt(std::span<const uint32_t> data, bool is_positive = true);

  /// @brief Construct from span of 16-bit words (big-endian order).
  /// @details Accepts any contiguous range: std::vector, std::array, C arrays, etc.
  explicit BigInt(std::span<const uint16_t> data, bool is_positive = true);

  /// @brief Construct from span of bytes (big-endian order).
  /// @details Accepts any contiguous range: std::vector, std::array, C arrays, etc.
  explicit BigInt(std::span<const uint8_t> data, bool is_positive = true);

  /// @brief Construct from vector of bits (big-endian order).
  /// @note std::vector<bool> is a special case that cannot use std::span.
  explicit BigInt(const std::vector<bool>& vec, bool is_positive = true);

  /// @brief Construct from any integral type.
  /// @details Handles int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t.
  template <std::integral T>
  explicit BigInt(T value);

  ~BigInt() = default;
  /// @}

  /// @name Assignment
  /// @{
  BigInt& operator=(const BigInt& other);
  BigInt& operator=(BigInt&& other) noexcept;
  /// @}

  /// @name Arithmetic Operators
  /// @{
  [[nodiscard]] BigInt operator+() const;                      ///< Unary plus
  [[nodiscard]] BigInt operator+(const BigInt& addend) const;  ///< Addition
  BigInt& operator+=(const BigInt& addend);                    ///< Addition assignment
  BigInt& operator++();                                        ///< Pre-increment
  BigInt operator++(int);                                      ///< Post-increment

  [[nodiscard]] BigInt operator-() const;                          ///< Unary minus (negation)
  [[nodiscard]] BigInt operator-(const BigInt& subtrahend) const;  ///< Subtraction
  BigInt& operator-=(const BigInt& subtrahend);                    ///< Subtraction assignment
  BigInt& operator--();                                            ///< Pre-decrement
  BigInt operator--(int);                                          ///< Post-decrement

  [[nodiscard]] BigInt operator*(const BigInt& multiplier) const;  ///< Multiplication
  BigInt& operator*=(const BigInt& multiplier);                    ///< Multiplication assignment

  /// @brief Integer division.
  /// @throws std::domain_error if divisor is zero.
  [[nodiscard]] BigInt operator/(const BigInt& divisor) const;
  /// @brief Division assignment.
  /// @throws std::domain_error if divisor is zero.
  BigInt& operator/=(const BigInt& divisor);

  /// @brief Modulo (remainder).
  /// @throws std::domain_error if divisor is zero.
  [[nodiscard]] BigInt operator%(const BigInt& divisor) const;

  /// @brief Fast modulo with small divisor using Horner's method.
  /// @details O(n) where n is number of 32-bit digits. No memory allocation.
  /// @param divisor The divisor (must be > 0).
  /// @return Remainder (always < divisor).
  /// @throws std::domain_error if divisor is zero.
  [[nodiscard]] uint32_t operator%(uint32_t divisor) const;

  /// @brief Modulo assignment.
  /// @throws std::domain_error if divisor is zero.
  BigInt& operator%=(const BigInt& divisor);
  /// @}

  /// @name Mathematical Functions
  /// @{

  /// @brief Compute base raised to exponent.
  /// @return base^exponent (0 if exponent is negative).
  friend BigInt pow(const BigInt& base, const BigInt& exponent);

  /// @brief Compute floor of log base 2.
  /// @return floor(log2(value)), or 0 if value <= 0.
  friend size_t log2(const BigInt& value) noexcept;

  /// @brief Compute (base^exponent) mod modulus efficiently.
  /// @throws std::domain_error if modulus is zero.
  friend BigInt powmod(const BigInt& base, const BigInt& exponent, const BigInt& modulus);

  /// @brief Compute modular multiplicative inverse.
  /// @return x such that (value * x) mod modulus == 1.
  /// @throws std::domain_error if modulus is zero or inverse doesn't exist.
  friend BigInt inversemod(BigInt value, const BigInt& modulus);

  /// @brief Check if a ≡ b (mod modulus).
  /// @throws std::domain_error if modulus is zero.
  friend bool congruencemod(const BigInt& a, const BigInt& b, const BigInt& modulus);

  /// @brief Check if two numbers are coprime (gcd == 1).
  friend bool isCoprime(const BigInt& a, const BigInt& b);

  /// @brief Compute Jacobi symbol (a/n).
  /// @return -1, 0, or 1 representing the Jacobi symbol.
  friend int8_t symbolJacobi(BigInt a, BigInt n);
  /// @}

  /// @name Bitwise Operators
  /// @{
  [[nodiscard]] BigInt operator~() const;                   ///< Bitwise NOT
  [[nodiscard]] BigInt operator&(const BigInt& rhs) const;  ///< Bitwise AND
  BigInt& operator&=(const BigInt& rhs);                    ///< Bitwise AND assignment
  [[nodiscard]] BigInt operator|(const BigInt& rhs) const;  ///< Bitwise OR
  BigInt& operator|=(const BigInt& rhs);                    ///< Bitwise OR assignment
  [[nodiscard]] BigInt operator^(const BigInt& rhs) const;  ///< Bitwise XOR
  BigInt& operator^=(const BigInt& rhs);                    ///< Bitwise XOR assignment

  [[nodiscard]] BigInt operator<<(size_t shift) const;  ///< Left shift
  BigInt& operator<<=(size_t shift);                    ///< Left shift assignment
  [[nodiscard]] BigInt operator>>(size_t shift) const;  ///< Right shift
  BigInt& operator>>=(size_t shift);                    ///< Right shift assignment

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

  /// @brief Compare with any integral type.
  /// @details Handles int8_t through int64_t and uint8_t through uint64_t.
  template <std::integral T>
  [[nodiscard]] std::strong_ordering operator<=>(T rhs) const noexcept;

  template <std::integral T>
  [[nodiscard]] bool operator==(T rhs) const noexcept;
  /// @}

  /// @name Utility Functions
  /// @{
  /// @brief Compute absolute value.
  friend BigInt abs(const BigInt& value);

  /// @brief Compute integer square root (floor).
  /// @throws std::domain_error if value is negative.
  friend BigInt sqrt(const BigInt& value);

  /// @brief Compute greatest common divisor.
  friend BigInt gcd(BigInt a, BigInt b);

  /// @brief Compute least common multiple.
  friend BigInt lcm(const BigInt& a, const BigInt& b);

  /// @brief Return the larger of two values.
  friend const BigInt& max(const BigInt& a, const BigInt& b) noexcept;

  /// @brief Return the smaller of two values.
  friend const BigInt& min(const BigInt& a, const BigInt& b) noexcept;

  /// @brief Swap two BigInt values.
  friend void swap(BigInt& lhs, BigInt& rhs) noexcept;
  /// @}

  /// @name Number Theory
  /// @{

  /// @brief Test if value is prime using Miller-Rabin.
  /// @details For numbers < 3,215,031,751, uses deterministic witnesses (100% accurate).
  ///          For larger numbers, uses probabilistic test with specified rounds.
  /// @param rounds Number of test rounds for large numbers (default 16).
  /// @return true if prime (or probably prime for large numbers), false if composite.
  [[nodiscard]] bool isProbablePrime(size_t rounds = 16) const;

  /// @brief Generate a random BigInt with specified number of bits.
  /// @param num_bits Number of bits in the result.
  /// @return Random BigInt in range [2^(num_bits-1), 2^num_bits - 1].
  [[nodiscard]] static BigInt randomBits(size_t num_bits);

  /// @brief Generate a random BigInt in range [0, upper_bound).
  /// @param upper_bound Upper bound (exclusive).
  /// @return Random BigInt in range [0, upper_bound).
  [[nodiscard]] static BigInt randomBelow(const BigInt& upper_bound);

  /// @brief Generate a random prime number with specified bit length.
  /// @details Uses randomBits() to generate candidates, filters by small primes,
  ///          then verifies with isProbablePrime(). Guaranteed to have exactly num_bits bits.
  /// @param num_bits Number of bits in the result (must be >= 2).
  /// @return Random prime with exactly num_bits bits.
  [[nodiscard]] static BigInt randomPrime(size_t num_bits);

  /// @brief Find the next prime >= this value.
  /// @return Smallest prime >= *this.
  [[nodiscard]] BigInt nextPrime() const;

  /// @}

  /// @name Stream I/O
  /// @{
  friend std::ostream& operator<<(std::ostream& out, const BigInt& value);
  friend std::istream& operator>>(std::istream& in, BigInt& value);
  /// @}

  /// @name Conversion
  /// @{

  /// @brief Convert to string in specified base (2, 10, or 16).
  /// @param base The numeric base.
  /// @note For formatted output with uppercase, showbase, etc., use operator<<
  ///       with standard stream manipulators (std::hex, std::uppercase, etc.)
  [[nodiscard]] std::string toStdString(uint8_t base = kDefaultOutputBase) const;

  /// @brief Convert to vector of bytes (big-endian order).
  /// @details Used for serialization, hashing, and crypto protocols.
  [[nodiscard]] explicit operator std::vector<uint8_t>() const;

  /// @brief Convert to any integral type (truncates if value doesn't fit).
  /// @details Handles int8_t through int64_t and uint8_t through uint64_t.
  ///          Signed types preserve sign, unsigned types return absolute value.
  template <std::integral T>
  [[nodiscard]] explicit operator T() const noexcept;

  [[nodiscard]] explicit operator bool() const noexcept;  ///< True if non-zero
  /// @}

  /// @name Queries
  /// @{
  [[nodiscard]] size_t bitLength() const noexcept;   ///< Number of bits (minimum 1 for zero)
  [[nodiscard]] size_t byteLength() const noexcept;  ///< Number of bytes needed
                                                     /// @}

 private:
  bool positive_;                 ///< Sign flag (true = positive or zero)
  std::vector<uint32_t> digits_;  ///< Little-endian 32-bit word storage

  /// @brief Compare magnitudes (absolute values) without allocation.
  /// @return -1 if |*this| < |other|, 0 if equal, 1 if |*this| > |other|
  [[nodiscard]] int compareMagnitude(const BigInt& other) const noexcept;

  BigInt operator*(uint32_t multiplier) const;
  BigInt& operator*=(uint32_t multiplier);
  std::pair<BigInt, BigInt> DivMod(const BigInt& divisor) const;
  friend BigInt BarrettReduction(const BigInt& dividend, const BigInt& divisor, const BigInt& mu);
  void alignTo(BigInt& other);
  void deleteZeroHighOrderDigit();
  BigInt shiftDigitsToHigh(size_t shift) const;
  BigInt shiftDigitsToLow(size_t shift) const;
  BigInt toBigIntDec() const;

  /// @brief Karatsuba multiplication helper for O(n^1.585) performance on large numbers.
  /// @param other The multiplier.
  /// @return Product of absolute values (sign handled by caller).
  [[nodiscard]] BigInt multiplyKaratsuba(const BigInt& other) const;

  /// @brief Schoolbook multiplication for small numbers or Karatsuba base case.
  /// @param other The multiplier.
  /// @return Product of absolute values (sign handled by caller).
  [[nodiscard]] BigInt multiplySchoolbook(const BigInt& other) const;

  // Helper for integral constructor/comparison
  void initFromUnsigned(uint64_t value, bool is_positive);
  void initFromSigned(int64_t value);
  [[nodiscard]] std::strong_ordering compareToSigned(int64_t rhs) const noexcept;
  [[nodiscard]] std::strong_ordering compareToUnsigned(uint64_t rhs) const noexcept;
  [[nodiscard]] uint64_t toUint64() const noexcept;
  [[nodiscard]] int64_t toInt64() const noexcept;
};

// Template implementations

template <std::integral T>
BigInt::BigInt(T value) {
  if constexpr (std::is_signed_v<T>) {
    initFromSigned(static_cast<int64_t>(value));
  } else {
    initFromUnsigned(static_cast<uint64_t>(value), true);
  }
}

template <std::integral T>
std::strong_ordering BigInt::operator<=>(T rhs) const noexcept {
  if constexpr (std::is_signed_v<T>) {
    return compareToSigned(static_cast<int64_t>(rhs));
  } else {
    return compareToUnsigned(static_cast<uint64_t>(rhs));
  }
}

template <std::integral T>
bool BigInt::operator==(T rhs) const noexcept {
  return (*this <=> rhs) == std::strong_ordering::equal;
}

template <std::integral T>
BigInt::operator T() const noexcept {
  if constexpr (std::is_signed_v<T>) {
    return static_cast<T>(toInt64());
  } else {
    return static_cast<T>(toUint64());
  }
}

}  // namespace bigint

namespace std {

/// @brief Hash specialization for BigInt to enable use in unordered containers.
template <>
struct hash<bigint::BigInt> {
  size_t operator()(const bigint::BigInt& value) const noexcept {
    // FNV-1a inspired hash combining sign and bytes
    size_t h = 14695981039346656037ULL;  // FNV offset basis
    h ^= static_cast<size_t>(value < 0);
    h *= 1099511628211ULL;  // FNV prime
    for (uint8_t byte : static_cast<std::vector<uint8_t>>(value)) {
      h ^= static_cast<size_t>(byte);
      h *= 1099511628211ULL;
    }
    return h;
  }
};

}  // namespace std

namespace bigint::literals {

/// @brief User-defined literal for BigInt.
/// @param str The string representation of the number.
/// @param len Length of the string.
/// @return BigInt constructed from the decimal string.
/// @example auto x = 123456789012345678901234567890_bigint;
inline BigInt operator""_bigint(const char* str, size_t len) {
  return BigInt(std::string(str, len), 10);
}

/// @brief User-defined literal for BigInt (integer form).
/// @param value The integer value.
/// @return BigInt constructed from the integer.
/// @example auto x = 12345_bigint;
inline BigInt operator""_bigint(unsigned long long value) {
  return BigInt(static_cast<uint64_t>(value));
}

}  // namespace bigint::literals
