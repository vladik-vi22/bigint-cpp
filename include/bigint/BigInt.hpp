/**
 * @file BigInt.hpp
 * @brief Arbitrary precision integer class.
 */
#pragma once

#include <bigint/Constants.hpp>

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

  /// @brief Assign from any integral type.
  template <std::integral T>
  BigInt& operator=(T value);
  /// @}

  /// @name Arithmetic Operators
  /// @{

  // Unary operators (members)
  [[nodiscard]] BigInt operator+() const;  ///< Unary plus
  [[nodiscard]] BigInt operator-() const;  ///< Unary minus (negation)

  // Increment/decrement (members)
  BigInt& operator++();    ///< Pre-increment
  BigInt operator++(int);  ///< Post-increment
  BigInt& operator--();    ///< Pre-decrement
  BigInt operator--(int);  ///< Post-decrement

  // Compound assignment (members)
  BigInt& operator+=(const BigInt& rhs);  ///< Addition assignment
  BigInt& operator-=(const BigInt& rhs);  ///< Subtraction assignment
  BigInt& operator*=(const BigInt& rhs);  ///< Multiplication assignment
  /// @throws std::domain_error if divisor is zero.
  BigInt& operator/=(const BigInt& rhs);  ///< Division assignment
  /// @throws std::domain_error if divisor is zero.
  BigInt& operator%=(const BigInt& rhs);  ///< Modulo assignment

  template <std::integral T>
  BigInt& operator+=(T rhs);
  template <std::integral T>
  BigInt& operator-=(T rhs);
  template <std::integral T>
  BigInt& operator*=(T rhs);
  template <std::integral T>
  BigInt& operator/=(T rhs);
  template <std::integral T>
  BigInt& operator%=(T rhs);

  // Binary arithmetic operators (friend functions for symmetry)
  /// @{
  friend BigInt operator+(const BigInt& lhs, const BigInt& rhs);
  friend BigInt operator-(const BigInt& lhs, const BigInt& rhs);
  friend BigInt operator*(const BigInt& lhs, const BigInt& rhs);
  /// @throws std::domain_error if divisor is zero.
  friend BigInt operator/(const BigInt& lhs, const BigInt& rhs);
  /// @throws std::domain_error if divisor is zero.
  friend BigInt operator%(const BigInt& lhs, const BigInt& rhs);

  /// @brief Arithmetic with integral types (BigInt op T).
  /// @details Optimized implementations for uint64_t/int64_t, smaller types cast up.
  /// @{
  template <std::integral T>
  friend BigInt operator+(const BigInt& lhs, T rhs);
  template <std::integral T>
  friend BigInt operator-(const BigInt& lhs, T rhs);
  template <std::integral T>
  friend BigInt operator*(const BigInt& lhs, T rhs);
  template <std::integral T>
  friend BigInt operator/(const BigInt& lhs, T rhs);
  template <std::integral T>
  friend T operator%(const BigInt& lhs, T rhs);
  /// @}

  /// @brief Reverse arithmetic operators (T op BigInt).
  /// @{
  template <std::integral T>
  friend BigInt operator+(T lhs, const BigInt& rhs) {
    return rhs + lhs;
  }
  template <std::integral T>
  friend BigInt operator-(T lhs, const BigInt& rhs) {
    return BigInt(lhs) - rhs;
  }
  template <std::integral T>
  friend BigInt operator*(T lhs, const BigInt& rhs) {
    return rhs * lhs;
  }
  template <std::integral T>
  friend BigInt operator/(T lhs, const BigInt& rhs) {
    return BigInt(lhs) / rhs;
  }
  template <std::integral T>
  friend T operator%(T lhs, const BigInt& rhs) {
    return static_cast<T>(BigInt(lhs) % rhs);
  }
  /// @}
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

  /// @name Logical Operators
  /// @{
  [[nodiscard]] bool operator!() const noexcept;                    ///< Logical NOT (true if zero)
  [[nodiscard]] bool operator&&(const BigInt& rhs) const noexcept;  ///< Logical AND
  [[nodiscard]] bool operator||(const BigInt& rhs) const noexcept;  ///< Logical OR
  /// @}

  /// @name Mathematical Functions
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

  /// @brief Compute base raised to exponent.
  /// @return base^exponent (0 if exponent is negative).
  friend BigInt pow(const BigInt& base, const BigInt& exponent);

  /// @brief Compute floor of log base 2.
  /// @return floor(log2(value)).
  /// @throws std::domain_error if value <= 0.
  friend size_t log2(const BigInt& value);

  /// @brief Compute quotient and remainder in one operation.
  /// @param dividend The dividend.
  /// @param divisor The divisor.
  /// @return Pair of (quotient, remainder).
  /// @throws std::domain_error if divisor is zero.
  friend std::pair<BigInt, BigInt> divmod(const BigInt& dividend, const BigInt& divisor);
  /// @}

  /// @name Modular Arithmetic
  /// @{

  /// @brief Compute (base^exponent) mod modulus efficiently.
  /// @details Uses Montgomery multiplication for large odd moduli,
  ///          Barrett reduction for medium moduli, standard method otherwise.
  /// @throws std::domain_error if modulus is zero.
  friend BigInt powmod(const BigInt& base, const BigInt& exponent, const BigInt& modulus);

  /// @brief Compute modular multiplicative inverse.
  /// @return x such that (value * x) mod modulus == 1.
  /// @throws std::domain_error if modulus is zero or inverse doesn't exist.
  friend BigInt inversemod(BigInt value, const BigInt& modulus);

  /// @brief Check if a ≡ b (mod modulus).
  /// @throws std::domain_error if modulus is zero.
  friend bool congruencemod(const BigInt& a, const BigInt& b, const BigInt& modulus);
  /// @}

  /// @name Number Theory
  /// @{

  /// @brief Check if two numbers are coprime (gcd == 1).
  friend bool isCoprime(const BigInt& a, const BigInt& b);

  /// @brief Compute Jacobi symbol (a/n).
  /// @return -1, 0, or 1 representing the Jacobi symbol.
  friend int8_t symbolJacobi(BigInt a, BigInt n);

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

  /// @name Bit Queries & Manipulation
  /// @{
  [[nodiscard]] size_t bitLength() const noexcept;      ///< Number of bits (minimum 1 for zero)
  [[nodiscard]] size_t byteLength() const noexcept;     ///< Number of bytes needed
  [[nodiscard]] bool testBit(size_t n) const noexcept;  ///< Test if bit at position n is set
  [[nodiscard]] size_t trailingZeros() const noexcept;  ///< Count trailing zero bits (0 for zero)
  [[nodiscard]] size_t popCount() const noexcept;       ///< Count number of set bits (1s)
  BigInt& setBit(size_t n);                             ///< Set bit at position n to 1
  BigInt& clearBit(size_t n);                           ///< Clear bit at position n to 0
  BigInt& flipBit(size_t n);                            ///< Toggle bit at position n
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

  /// @name Stream I/O
  /// @{
  friend std::ostream& operator<<(std::ostream& out, const BigInt& value);
  friend std::istream& operator>>(std::istream& in, BigInt& value);
  /// @}

  /// @name Utility
  /// @{

  /// @brief Return the larger of two values.
  friend const BigInt& max(const BigInt& a, const BigInt& b) noexcept;

  /// @brief Return the smaller of two values.
  friend const BigInt& min(const BigInt& a, const BigInt& b) noexcept;

  /// @brief Swap two BigInt values.
  friend void swap(BigInt& lhs, BigInt& rhs) noexcept;
  /// @}

 private:
  bool positive_;                 ///< Sign flag (true = positive or zero)
  std::vector<uint32_t> digits_;  ///< Little-endian 32-bit word storage

  /// @brief Compare magnitudes (absolute values) without allocation.
  /// @return -1 if |*this| < |other|, 0 if equal, 1 if |*this| > |other|
  [[nodiscard]] int compareMagnitude(const BigInt& other) const noexcept;

  [[nodiscard]] std::pair<BigInt, BigInt> divmodSingleWord(uint32_t divisor,
                                                           bool quotient_positive) const;
  [[nodiscard]] std::pair<BigInt, BigInt> divmodKnuth(const BigInt& divisor, size_t m, size_t n,
                                                      bool quotient_positive) const;
  void deleteZeroHighOrderDigit();
  [[nodiscard]] BigInt shiftDigitsToHigh(size_t shift) const;
  [[nodiscard]] BigInt shiftDigitsToLow(size_t shift) const;
  [[nodiscard]] BigInt toBigIntDec() const;

  /// @brief Toom-Cook 3-way multiplication for O(n^1.465) performance on large numbers.
  /// @param other The multiplier.
  /// @return Product of absolute values (sign handled by caller).
  /// @details Splits operands into 3 parts, evaluates at 5 points (0, 1, -1, 2, ∞),
  ///          performs 5 recursive multiplications, then interpolates the result.
  [[nodiscard]] BigInt multiplyToom3(const BigInt& other) const;

  /// @brief Karatsuba multiplication helper for O(n^1.585) performance on medium numbers.
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

  // Integral arithmetic implementations (called by templates)
  [[nodiscard]] BigInt addInt64(int64_t rhs) const;
  [[nodiscard]] BigInt addUint64(uint64_t rhs) const;
  [[nodiscard]] BigInt subInt64(int64_t rhs) const;
  [[nodiscard]] BigInt subUint64(uint64_t rhs) const;
  [[nodiscard]] BigInt mulInt64(int64_t rhs) const;
  [[nodiscard]] BigInt mulUint64(uint64_t rhs) const;
  [[nodiscard]] BigInt mulUint32(uint32_t rhs) const;  // used by mulUint64
  [[nodiscard]] BigInt divInt64(int64_t rhs) const;
  [[nodiscard]] BigInt divUint64(uint64_t rhs) const;
  [[nodiscard]] int64_t modInt64(int64_t rhs) const;
  [[nodiscard]] uint64_t modUint64(uint64_t rhs) const;
  [[nodiscard]] uint32_t modUint32(uint32_t rhs) const;  // used by modUint64
};

// Template implementations

template <std::integral T>
BigInt::BigInt(T value) : positive_(true), digits_() {
  if constexpr (std::is_signed_v<T>) {
    initFromSigned(static_cast<int64_t>(value));
  } else {
    initFromUnsigned(static_cast<uint64_t>(value), true);
  }
}

template <std::integral T>
BigInt& BigInt::operator=(T value) {
  *this = BigInt(value);
  return *this;
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

template <std::integral T>
BigInt operator+(const BigInt& lhs, T rhs) {
  if constexpr (std::is_signed_v<T>) {
    return lhs.addInt64(static_cast<int64_t>(rhs));
  } else {
    return lhs.addUint64(static_cast<uint64_t>(rhs));
  }
}

template <std::integral T>
BigInt operator-(const BigInt& lhs, T rhs) {
  if constexpr (std::is_signed_v<T>) {
    return lhs.subInt64(static_cast<int64_t>(rhs));
  } else {
    return lhs.subUint64(static_cast<uint64_t>(rhs));
  }
}

template <std::integral T>
BigInt operator*(const BigInt& lhs, T rhs) {
  if constexpr (std::is_signed_v<T>) {
    return lhs.mulInt64(static_cast<int64_t>(rhs));
  } else {
    return lhs.mulUint64(static_cast<uint64_t>(rhs));
  }
}

template <std::integral T>
BigInt operator/(const BigInt& lhs, T rhs) {
  if constexpr (std::is_signed_v<T>) {
    return lhs.divInt64(static_cast<int64_t>(rhs));
  } else {
    return lhs.divUint64(static_cast<uint64_t>(rhs));
  }
}

template <std::integral T>
T operator%(const BigInt& lhs, T rhs) {
  if constexpr (std::is_signed_v<T>) {
    return static_cast<T>(lhs.modInt64(static_cast<int64_t>(rhs)));
  } else {
    return static_cast<T>(lhs.modUint64(static_cast<uint64_t>(rhs)));
  }
}

template <std::integral T>
BigInt& BigInt::operator+=(T rhs) {
  *this = *this + rhs;
  return *this;
}

template <std::integral T>
BigInt& BigInt::operator-=(T rhs) {
  *this = *this - rhs;
  return *this;
}

template <std::integral T>
BigInt& BigInt::operator*=(T rhs) {
  *this = *this * rhs;
  return *this;
}

template <std::integral T>
BigInt& BigInt::operator/=(T rhs) {
  *this = *this / rhs;
  return *this;
}

template <std::integral T>
BigInt& BigInt::operator%=(T rhs) {
  *this = *this % BigInt(rhs);
  return *this;
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
