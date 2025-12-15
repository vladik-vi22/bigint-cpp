/**
 * @file BigInt.cpp
 * @brief Implementation of arbitrary precision integer arithmetic.
 *
 * @details
 * ## Internal Representation
 * - **Storage**: Little-endian vector of 32-bit words (`digits_`)
 *   - `digits_[0]` contains the least significant 32 bits
 *   - `digits_.back()` contains the most significant bits
 * - **Sign**: Sign-magnitude encoding via `positive_` flag
 *   - Zero is always represented as positive
 *
 * ## Algorithms
 * - **Multiplication**: Schoolbook for small numbers, Karatsuba for large (threshold: 32 words)
 * - **Division**: Knuth Algorithm D (multi-word long division with normalization)
 * - **Modular Exponentiation**: Montgomery CIOS for large odd moduli, standard square-and-multiply
 * otherwise
 * - **Primality**: Miller-Rabin with deterministic witnesses for small numbers
 * - **GCD**: Euclidean algorithm with fast division (Knuth Algorithm D)
 *
 * ## Performance Considerations
 * - Operations reserve capacity to minimize reallocations
 * - Karatsuba multiplication reduces complexity from O(n²) to O(n^1.585)
 * - Montgomery multiplication avoids expensive division in modular exponentiation
 * - Small prime sieve accelerates primality testing
 */

#include <bigint/BigInt.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cmath>
#include <compare>
#include <iomanip>
#include <iterator>
#include <random>
#include <ranges>
#include <sstream>
#include <stdexcept>

#include "BarrettContext.hpp"
#include "MontgomeryContext.hpp"

namespace bigint {

// ============================================================================
// Internal Constants
// ============================================================================

namespace {

/// Number of bits per digit (32-bit words)
constexpr size_t kBitsPerDigit = 32;

/// Mask for extracting bit position within a digit (0-31)
constexpr size_t kBitIndexMask = kBitsPerDigit - 1;

/// Number of bits to shift to convert bit index to digit index
constexpr size_t kDigitShift = 5;  // log2(32)

/// Base for carry calculations (2^32)
constexpr uint64_t kDigitBase = 1ULL << kBitsPerDigit;

/// Base for decimal conversion (10^9, largest power of 10 fitting in 32 bits)
constexpr uint32_t kDecimalBase = 1'000'000'000;

/// Digits per decimal cell (log10(kDecimalBase))
constexpr uint8_t kDecimalCellSize = 9;

/// Threshold for switching from schoolbook to Karatsuba multiplication
constexpr size_t kKaratsubaThreshold = 32;

/// Threshold for switching from Karatsuba to Toom-Cook 3-way multiplication
/// Higher threshold due to interpolation overhead (division by 3, many additions)
constexpr size_t kToom3Threshold = 256;

/// Small primes for quick divisibility rejection in prime generation
constexpr auto kSmallPrimes =
    std::to_array<uint32_t>({3,  5,  7,  11, 13, 17, 19, 23, 29, 31,  37,  41,  43,  47, 53,
                             59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113});

/**
 * @brief Converts a decimal string to a binary string.
 * @param decimal_str Decimal number as string (no sign prefix).
 * @return Binary representation as string.
 * @details Uses repeated division by 2, processing 9 decimal digits at a time
 *          for efficiency. Each iteration extracts one binary digit.
 */
std::string decimalToBinaryString(std::string decimal_str) {
  if (decimal_str == "0") {
    return "0";
  }

  std::string binary_str;
  std::vector<uint32_t> cells;
  uint32_t carry_next = 0;
  uint32_t carry_current = 0;

  // Pad to multiple of kDecimalCellSize for uniform cell processing
  while (decimal_str.length() % kDecimalCellSize) {
    decimal_str.insert(0, 1, '0');
  }

  // Parse decimal string into cells of 9 digits each
  const size_t num_cells = decimal_str.length() / kDecimalCellSize;
  cells.reserve(num_cells);
  for (size_t i = 0; i < num_cells; ++i) {
    cells.emplace_back(static_cast<uint32_t>(
        std::stoul(decimal_str.substr(i * kDecimalCellSize, kDecimalCellSize), nullptr, 10)));
  }

  // Repeated division by 2 to extract binary digits (build in reverse, then flip)
  while (!std::ranges::all_of(cells, [](uint32_t c) { return c == 0; })) {
    carry_next = 0;
    for (auto& cell : cells) {
      carry_current = carry_next;
      carry_next = (cell & 1);
      cell = (cell + carry_current * kDecimalBase) >> 1;
    }
    binary_str.push_back(carry_next ? '1' : '0');
  }

  std::ranges::reverse(binary_str);
  return binary_str;
}

// ============================================================================
// Division Algorithm Helpers (Knuth's Algorithm D)
// ============================================================================

/// Word and double-word types for division
using Word = uint32_t;
using DWord = uint64_t;
using SignedDWord = int64_t;

/// Mask for extracting lower 32 bits
constexpr DWord kWordMask = 0xFFFFFFFF;

/**
 * @brief Counts leading zeros in a 32-bit word.
 * @param word The word to analyze.
 * @return Number of leading zero bits (0-32).
 */
[[nodiscard]] inline size_t countLeadingZeros(Word word) noexcept {
  return static_cast<size_t>(std::countl_zero(word));
}

/**
 * @brief Shifts a word vector left by a given number of bits.
 * @param src Source vector.
 * @param dst Destination vector (must be pre-sized).
 * @param shift Number of bits to shift (0-31).
 */
inline void shiftWordsLeft(const std::vector<Word>& src, std::vector<Word>& dst, size_t shift) {
  if (shift == 0) {
    std::ranges::copy(src, dst.begin());
    return;
  }
  Word carry = 0;
  const size_t right_shift = kBitsPerDigit - shift;
  for (size_t i = 0; i < src.size(); ++i) {
    dst[i] = (src[i] << shift) | carry;
    carry = src[i] >> right_shift;
  }
  if (dst.size() > src.size()) {
    dst[src.size()] = carry;
  }
}

/**
 * @brief Shifts a word vector right by a given number of bits.
 * @param src Source vector.
 * @param dst Destination vector (must be pre-sized).
 * @param n Number of words to process.
 * @param shift Number of bits to shift (0-31).
 */
inline void shiftWordsRight(const std::vector<Word>& src, std::vector<Word>& dst, size_t n,
                            size_t shift) {
  if (shift == 0) {
    std::ranges::copy_n(src.begin(), n, dst.begin());
    return;
  }
  Word carry = 0;
  const size_t left_shift = kBitsPerDigit - shift;
  for (size_t i = n; i > 0; --i) {
    dst[i - 1] = (src[i - 1] >> shift) | carry;
    carry = src[i - 1] << left_shift;
  }
}

/**
 * @brief Computes trial quotient digit using Knuth's method.
 * @param u_high High word of dividend segment.
 * @param u_mid Middle word of dividend segment.
 * @param u_low Low word of dividend segment.
 * @param v_high High word of divisor.
 * @param v_mid Second-highest word of divisor.
 * @return Refined trial quotient (may still be 1 too large).
 */
[[nodiscard]] inline DWord computeTrialQuotient(Word u_high, Word u_mid, Word u_low, DWord v_high,
                                                DWord v_mid) noexcept {
  // Form two-word dividend
  const DWord u_combined = (static_cast<DWord>(u_high) << kBitsPerDigit) | u_mid;
  DWord q_hat = u_combined / v_high;
  DWord r_hat = u_combined % v_high;

  // Refine using second digit of divisor
  while (q_hat >= kDigitBase || q_hat * v_mid > (r_hat << kBitsPerDigit) + u_low) {
    --q_hat;
    r_hat += v_high;
    if (r_hat >= kDigitBase) {
      break;
    }
  }
  return q_hat;
}

/**
 * @brief Multiplies divisor by q_hat and subtracts from dividend segment.
 * @param u Dividend words (modified in place).
 * @param v Divisor words.
 * @param q_hat Trial quotient digit.
 * @param idx Starting index in u.
 * @param n Number of divisor words.
 * @return Final borrow (negative if subtraction underflowed).
 */
inline SignedDWord multiplyAndSubtract(std::vector<Word>& u, const std::vector<Word>& v,
                                       DWord q_hat, size_t idx, size_t n) {
  SignedDWord borrow = 0;
  for (size_t i = 0; i < n; ++i) {
    const DWord product = q_hat * v[i];
    const SignedDWord diff = static_cast<SignedDWord>(u[idx + i]) -
                             static_cast<SignedDWord>(product & kWordMask) + borrow;
    u[idx + i] = static_cast<Word>(diff & kWordMask);
    borrow = (diff >> kBitsPerDigit) - static_cast<SignedDWord>(product >> kBitsPerDigit);
  }
  const SignedDWord final_diff = static_cast<SignedDWord>(u[idx + n]) + borrow;
  u[idx + n] = static_cast<Word>(final_diff & kWordMask);
  return final_diff;
}

/**
 * @brief Adds divisor back to dividend segment (rare correction step).
 * @param u Dividend words (modified in place).
 * @param v Divisor words.
 * @param idx Starting index in u.
 * @param n Number of divisor words.
 * @details Called when trial quotient was 1 too large (probability ~2/base).
 */
inline void addBack(std::vector<Word>& u, const std::vector<Word>& v, size_t idx, size_t n) {
  DWord carry = 0;
  for (size_t i = 0; i < n; ++i) {
    const DWord sum = static_cast<DWord>(u[idx + i]) + v[i] + carry;
    u[idx + i] = static_cast<Word>(sum & kWordMask);
    carry = sum >> kBitsPerDigit;
  }
  u[idx + n] += static_cast<Word>(carry);
}

}  // anonymous namespace

// ============================================================================
// Constructors
// ============================================================================

BigInt::BigInt() : positive_(true), digits_() {}

BigInt::BigInt(const BigInt& other) : positive_(other.positive_), digits_(other.digits_) {}

BigInt::BigInt(BigInt&& other) noexcept
    : positive_(other.positive_), digits_(std::move(other.digits_)) {
  other.positive_ = true;
}

BigInt::BigInt(std::string str, const uint8_t base) {
  // Handle empty string - treat as zero
  if (str.empty()) {
    positive_ = true;
    digits_.push_back(0);
    return;
  }

  const uint8_t cell_size =
      base == kBaseHexadecimal ? (sizeof(uint32_t) * 2) : (sizeof(uint32_t) * 8);
  if (str[0] == '-') {
    positive_ = false;
    str.erase(0, 1);
    // Handle "-" alone as zero
    if (str.empty()) {
      positive_ = true;
      digits_.push_back(0);
      return;
    }
  } else {
    positive_ = true;
  }

  // Validate input characters before processing
  auto isValidChar = [base](char c) -> bool {
    if (base == kBaseBinary) {
      return c == '0' || c == '1';
    } else if (base == kBaseDecimal) {
      return c >= '0' && c <= '9';
    } else if (base == kBaseHexadecimal) {
      return (c >= '0' && c <= '9') || (c >= 'A' && c <= 'F') || (c >= 'a' && c <= 'f');
    }
    return false;
  };

  for (size_t i = 0; i < str.size(); ++i) {
    if (!isValidChar(str[i])) {
      throw std::invalid_argument("Invalid character '" + std::string(1, str[i]) +
                                  "' at position " + std::to_string(i) + " for base " +
                                  std::to_string(base));
    }
  }

  if (base == kBaseDecimal) {
    str = decimalToBinaryString(str);
  }
  while (str.length() % cell_size) {
    str.insert(0, 1, '0');
  }
  size_t num_cells = str.length() / cell_size;
  digits_.reserve(num_cells);
  for (size_t i = 0; i < num_cells; ++i) {
    digits_.emplace(digits_.begin(),
                    static_cast<uint32_t>(
                        std::stoul(str.substr(i * cell_size, cell_size), nullptr,
                                   base == kBaseHexadecimal ? kBaseHexadecimal : kBaseBinary)));
  }
}

BigInt::BigInt(const std::vector<bool>& vec, const bool is_positive_) {
  const size_t remainder_bits = vec.size() & kBitIndexMask;
  const size_t full_digits = vec.size() >> kDigitShift;
  digits_.reserve(remainder_bits ? full_digits + 1 : full_digits);

  uint32_t digit = 0;
  auto it = vec.crbegin();

  // Process full 32-bit digits
  for (size_t i = 0; i < full_digits; ++i) {
    digit = 0;
    for (size_t bit_idx = 0; bit_idx < kBitsPerDigit; ++bit_idx) {
      digit |= static_cast<uint32_t>(*it) << bit_idx;
      ++it;
    }
    digits_.emplace_back(digit);
  }

  // Process remaining bits (partial digit)
  if (remainder_bits) {
    digit = 0;
    for (size_t bit_idx = 0; bit_idx < remainder_bits; ++bit_idx) {
      digit |= static_cast<uint32_t>(*it) << bit_idx;
      ++it;
    }
    digits_.emplace_back(digit);
  }

  deleteZeroHighOrderDigit();
  positive_ = is_positive_;
}

BigInt::BigInt(std::span<const uint32_t> data, const bool is_positive_) : positive_(is_positive_) {
  digits_.reserve(data.size());
  // Input is big-endian, internal storage is little-endian
  std::ranges::copy(std::views::reverse(data), std::back_inserter(digits_));
  deleteZeroHighOrderDigit();
}

BigInt::BigInt(std::span<const uint16_t> data, const bool is_positive_) {
  digits_.reserve((data.size() + 1) / 2);
  // Process from end (big-endian input) in pairs
  size_t i = data.size();
  while (i >= 2) {
    i -= 2;
    digits_.emplace_back(static_cast<uint32_t>(data[i + 1]) |
                         (static_cast<uint32_t>(data[i]) << 16));
  }
  if (i == 1) {
    digits_.emplace_back(static_cast<uint32_t>(data[0]));
  }
  deleteZeroHighOrderDigit();
  positive_ = is_positive_;
}

BigInt::BigInt(std::span<const uint8_t> data, const bool is_positive_) {
  digits_.reserve((data.size() + 3) / 4);
  // Process from end (big-endian input) in groups of 4
  size_t i = data.size();
  while (i >= 4) {
    i -= 4;
    digits_.emplace_back(
        static_cast<uint32_t>(data[i + 3]) | (static_cast<uint32_t>(data[i + 2]) << 8) |
        (static_cast<uint32_t>(data[i + 1]) << 16) | (static_cast<uint32_t>(data[i]) << 24));
  }
  // Handle remaining bytes
  if (i == 3) {
    digits_.emplace_back(static_cast<uint32_t>(data[2]) | (static_cast<uint32_t>(data[1]) << 8) |
                         (static_cast<uint32_t>(data[0]) << 16));
  } else if (i == 2) {
    digits_.emplace_back(static_cast<uint32_t>(data[1]) | (static_cast<uint32_t>(data[0]) << 8));
  } else if (i == 1) {
    digits_.emplace_back(static_cast<uint32_t>(data[0]));
  }
  deleteZeroHighOrderDigit();
  positive_ = is_positive_;
}

void BigInt::initFromUnsigned(const uint64_t value, const bool is_positive) {
  positive_ = is_positive;
  digits_.reserve(2);
  digits_.emplace_back(static_cast<uint32_t>(value));
  digits_.emplace_back(static_cast<uint32_t>(value >> kBitsPerDigit));
  deleteZeroHighOrderDigit();
}

void BigInt::initFromSigned(const int64_t value) {
  positive_ = value >= 0;
  const auto abs_value = static_cast<uint64_t>(value >= 0 ? value : -value);
  digits_.reserve(2);
  digits_.emplace_back(static_cast<uint32_t>(abs_value));
  digits_.emplace_back(static_cast<uint32_t>(abs_value >> kBitsPerDigit));
  deleteZeroHighOrderDigit();
}

// ============================================================================
// Assignment Operators
// ============================================================================

BigInt& BigInt::operator=(const BigInt& other) {
  if (this != &other) {
    digits_ = other.digits_;
    positive_ = other.positive_;
  }
  return *this;
}

BigInt& BigInt::operator=(BigInt&& other) noexcept {
  if (this != &other) {
    digits_ = std::move(other.digits_);
    positive_ = other.positive_;
    other.positive_ = true;
  }
  return *this;
}

// ============================================================================
// Arithmetic Operators
// ============================================================================

BigInt BigInt::operator+() const {
  return *this;
}

BigInt BigInt::operator+(const BigInt& addend) const {
  if (positive_ && addend.positive_) {
    BigInt sum;
    uint32_t carry = 0;
    uint64_t temp_sum;
    const bool this_larger = (digits_.size() >= addend.digits_.size());

    sum.digits_.reserve(this_larger ? digits_.size() + 1 : addend.digits_.size() + 1);

    auto larger_it = this_larger ? digits_.cbegin() : addend.digits_.cbegin();
    auto smaller_it = this_larger ? addend.digits_.cbegin() : digits_.cbegin();
    auto larger_end = this_larger ? digits_.cend() : addend.digits_.cend();
    auto smaller_end = this_larger ? addend.digits_.cend() : digits_.cend();

    while (smaller_it != smaller_end) {
      temp_sum = static_cast<uint64_t>(*larger_it) + static_cast<uint64_t>(*smaller_it) +
                 static_cast<uint64_t>(carry);
      sum.digits_.emplace_back(static_cast<uint32_t>(temp_sum));
      carry = static_cast<uint32_t>(temp_sum >> kBitsPerDigit);
      ++larger_it;
      ++smaller_it;
    }

    while (larger_it != larger_end) {
      temp_sum = static_cast<uint64_t>(*larger_it) + static_cast<uint64_t>(carry);
      sum.digits_.emplace_back(static_cast<uint32_t>(temp_sum));
      carry = static_cast<uint32_t>(temp_sum >> kBitsPerDigit);
      ++larger_it;
    }

    if (carry) {
      sum.digits_.emplace_back(carry);
    }
    sum.positive_ = true;
    return sum;
  } else if (positive_ && !addend.positive_) {
    return *this - abs(addend);
  } else if (!positive_ && addend.positive_) {
    return addend - abs(*this);
  } else {  // !positive_ && !addend.positive_
    return -(abs(*this) + abs(addend));
  }
}

BigInt& BigInt::operator+=(const BigInt& addend) {
  // Handle zero cases
  if (addend == 0) {
    return *this;
  }
  if (*this == 0) {
    *this = addend;
    return *this;
  }

  // Same sign: add magnitudes
  if (positive_ == addend.positive_) {
    if (digits_.size() < addend.digits_.size()) {
      digits_.resize(addend.digits_.size(), 0);
    }

    uint32_t carry = 0;
    size_t i = 0;
    for (; i < addend.digits_.size(); ++i) {
      uint64_t sum = static_cast<uint64_t>(digits_[i]) + static_cast<uint64_t>(addend.digits_[i]) +
                     static_cast<uint64_t>(carry);
      digits_[i] = static_cast<uint32_t>(sum);
      carry = static_cast<uint32_t>(sum >> kBitsPerDigit);
    }
    // Propagate carry through remaining digits
    for (; carry && i < digits_.size(); ++i) {
      uint64_t sum = static_cast<uint64_t>(digits_[i]) + carry;
      digits_[i] = static_cast<uint32_t>(sum);
      carry = static_cast<uint32_t>(sum >> kBitsPerDigit);
    }
    if (carry) {
      digits_.push_back(carry);
    }
  } else {
    // Different signs: subtract smaller magnitude from larger
    int mag_cmp = compareMagnitude(addend);
    if (mag_cmp == 0) {
      digits_.clear();
      digits_.push_back(0);
      positive_ = true;
      return *this;
    }

    const BigInt* larger = (mag_cmp > 0) ? this : &addend;
    const BigInt* smaller = (mag_cmp > 0) ? &addend : this;
    bool result_positive = (mag_cmp > 0) ? positive_ : addend.positive_;

    std::vector<uint32_t> result(larger->digits_.size(), 0);
    uint8_t borrow = 0;

    for (size_t i = 0; i < larger->digits_.size(); ++i) {
      int64_t smaller_digit =
          (i < smaller->digits_.size()) ? static_cast<int64_t>(smaller->digits_[i]) : 0;
      int64_t diff = static_cast<int64_t>(larger->digits_[i]) - smaller_digit - borrow;
      if (diff >= 0) {
        result[i] = static_cast<uint32_t>(diff);
        borrow = 0;
      } else {
        result[i] = static_cast<uint32_t>(diff + static_cast<int64_t>(kDigitBase));
        borrow = 1;
      }
    }

    digits_ = std::move(result);
    positive_ = result_positive;
    deleteZeroHighOrderDigit();
  }

  return *this;
}

BigInt& BigInt::operator++() {
  if (positive_) {
    if (digits_.empty()) {
      digits_.push_back(1);
      return *this;
    }
    for (size_t i = 0; i < digits_.size(); ++i) {
      if (digits_[i] < UINT32_MAX) {
        ++digits_[i];
        return *this;
      }
      digits_[i] = 0;
    }
    digits_.push_back(1);
  } else {
    for (size_t i = 0; i < digits_.size(); ++i) {
      if (digits_[i] > 0) {
        --digits_[i];
        deleteZeroHighOrderDigit();
        if (*this == 0) {
          positive_ = true;
        }
        return *this;
      }
      digits_[i] = UINT32_MAX;
    }
  }
  return *this;
}

BigInt BigInt::operator++(int) {
  const BigInt old_value = *this;
  ++(*this);
  return old_value;
}

BigInt BigInt::operator-() const {
  BigInt negated = *this;
  negated.positive_ = !positive_;
  return negated;
}

BigInt BigInt::operator-(const BigInt& subtrahend) const {
  if (positive_ && subtrahend.positive_) {
    if (*this >= subtrahend) {
      BigInt diff;
      uint8_t borrow = 0;
      diff.digits_.reserve(digits_.size());

      auto minuend_it = digits_.cbegin();
      auto subtrahend_it = subtrahend.digits_.cbegin();

      while (subtrahend_it != subtrahend.digits_.cend()) {
        int64_t temp_diff = static_cast<int64_t>(*minuend_it) -
                            static_cast<int64_t>(*subtrahend_it) - static_cast<int64_t>(borrow);
        if (temp_diff >= 0) {
          diff.digits_.emplace_back(static_cast<uint32_t>(temp_diff));
          borrow = 0;
        } else {
          diff.digits_.emplace_back(
              static_cast<uint32_t>(temp_diff + static_cast<int64_t>(kDigitBase)));
          borrow = 1;
        }
        ++minuend_it;
        ++subtrahend_it;
      }

      while (minuend_it != digits_.cend()) {
        int64_t temp_diff = static_cast<int64_t>(*minuend_it) - static_cast<int64_t>(borrow);
        if (temp_diff >= 0) {
          diff.digits_.emplace_back(static_cast<uint32_t>(temp_diff));
          borrow = 0;
        } else {
          diff.digits_.emplace_back(
              static_cast<uint32_t>(temp_diff + static_cast<int64_t>(kDigitBase)));
          borrow = 1;
        }
        ++minuend_it;
      }

      diff.deleteZeroHighOrderDigit();
      diff.positive_ = true;
      return diff;
    } else {
      return -(subtrahend - *this);
    }
  } else if (!positive_ && subtrahend.positive_) {
    return -(abs(*this) + subtrahend);
  } else if (positive_ && !subtrahend.positive_) {
    return *this + abs(subtrahend);
  } else {  // !positive_ && !subtrahend.positive_
    return abs(subtrahend) - abs(*this);
  }
}

BigInt& BigInt::operator-=(const BigInt& subtrahend) {
  if (subtrahend == 0) {
    return *this;
  }
  BigInt negated = subtrahend;
  negated.positive_ = !subtrahend.positive_;
  return *this += negated;
}

BigInt& BigInt::operator--() {
  if (positive_) {
    if (digits_.empty() || *this == 0) {
      digits_.clear();
      digits_.push_back(1);
      positive_ = false;
      return *this;
    }
    for (size_t i = 0; i < digits_.size(); ++i) {
      if (digits_[i] > 0) {
        --digits_[i];
        deleteZeroHighOrderDigit();
        return *this;
      }
      digits_[i] = UINT32_MAX;
    }
  } else {
    if (digits_.empty()) {
      digits_.push_back(1);
      return *this;
    }
    for (size_t i = 0; i < digits_.size(); ++i) {
      if (digits_[i] < UINT32_MAX) {
        ++digits_[i];
        return *this;
      }
      digits_[i] = 0;
    }
    digits_.push_back(1);
  }
  return *this;
}

BigInt BigInt::operator--(int) {
  const BigInt old_value = *this;
  --(*this);
  return old_value;
}

BigInt BigInt::operator*(uint32_t multiplier) const {
  BigInt product;
  uint32_t carry = 0;
  product.digits_.reserve(digits_.size() + 1);

  for (const auto& digit : digits_) {
    uint64_t temp_product =
        static_cast<uint64_t>(digit) * static_cast<uint64_t>(multiplier) + carry;
    product.digits_.emplace_back(static_cast<uint32_t>(temp_product));
    carry = static_cast<uint32_t>(temp_product >> kBitsPerDigit);
  }

  if (carry) {
    product.digits_.emplace_back(carry);
  }
  product.positive_ = positive_;
  return product;
}

BigInt& BigInt::operator*=(uint32_t multiplier) {
  uint32_t carry = 0;
  for (auto& digit : digits_) {
    uint64_t product = static_cast<uint64_t>(digit) * multiplier + carry;
    digit = static_cast<uint32_t>(product);
    carry = static_cast<uint32_t>(product >> kBitsPerDigit);
  }
  if (carry) {
    digits_.emplace_back(carry);
  }
  return *this;
}

BigInt BigInt::operator*(const BigInt& multiplier) const {
  if (*this == 0 || multiplier == 0) {
    return BigInt();
  }

  BigInt product;
  const size_t max_size = std::max(digits_.size(), multiplier.digits_.size());

  if (max_size >= kToom3Threshold) {
    product = multiplyToom3(multiplier);
  } else if (max_size >= kKaratsubaThreshold) {
    product = multiplyKaratsuba(multiplier);
  } else {
    product = multiplySchoolbook(multiplier);
  }

  product.positive_ = (positive_ == multiplier.positive_);
  return product;
}

BigInt& BigInt::operator*=(const BigInt& multiplier) {
  *this = *this * multiplier;
  return *this;
}

std::pair<BigInt, BigInt> BigInt::divmod(const BigInt& divisor) const {
  if (divisor == 0) {
    throw std::domain_error("Division by zero");
  }

  const bool quotient_positive = (positive_ == divisor.positive_);

  // Compare magnitudes to handle trivial cases
  const int cmp = compareMagnitude(divisor);
  if (cmp < 0) {
    // |dividend| < |divisor| → quotient = 0, remainder = dividend
    BigInt remainder = *this;
    remainder.positive_ = positive_ || (remainder == 0);
    return {BigInt(), remainder};
  }
  if (cmp == 0) {
    // |dividend| == |divisor| → quotient = ±1, remainder = 0
    BigInt quotient(1);
    quotient.positive_ = quotient_positive;
    return {quotient, BigInt()};
  }

  const size_t m = digits_.size();
  const size_t n = divisor.digits_.size();

  // Fast path: single-word divisor
  if (n == 1) {
    return divmodSingleWord(divisor.digits_[0], quotient_positive);
  }

  // Knuth's Algorithm D (TAOCP Vol 2, Section 4.3.1)
  return divmodKnuth(divisor, m, n, quotient_positive);
}

/**
 * @brief Fast division by single-word divisor.
 */
std::pair<BigInt, BigInt> BigInt::divmodSingleWord(uint32_t d, bool quotient_positive) const {
  const size_t m = digits_.size();
  BigInt quotient;
  quotient.digits_.resize(m);

  uint64_t carry = 0;
  for (size_t i = m; i-- > 0;) {
    const uint64_t cur = carry * kDigitBase + digits_[i];
    quotient.digits_[i] = static_cast<uint32_t>(cur / d);
    carry = cur % d;
  }

  quotient.deleteZeroHighOrderDigit();
  quotient.positive_ = quotient_positive || !quotient;

  BigInt remainder(static_cast<uint32_t>(carry));
  remainder.positive_ = positive_ || !remainder;
  return {std::move(quotient), std::move(remainder)};
}

/**
 * @brief Multi-word division using Knuth's Algorithm D.
 */
std::pair<BigInt, BigInt> BigInt::divmodKnuth(const BigInt& divisor, size_t m, size_t n,
                                              bool quotient_positive) const {
  // D1: Normalize - shift so that v[n-1] >= base/2
  const size_t shift = countLeadingZeros(divisor.digits_[n - 1]);

  // Create normalized copies: u has m+1 words, v has n words
  std::vector<Word> u(m + 1, 0);
  std::vector<Word> v(n, 0);
  shiftWordsLeft(divisor.digits_, v, shift);
  shiftWordsLeft(digits_, u, shift);

  // D2-D7: Main loop - compute quotient digits from high to low
  const size_t q_size = m - n + 1;
  std::vector<Word> q(q_size, 0);

  const DWord v_high = v[n - 1];
  const DWord v_mid = (n >= 2) ? v[n - 2] : 0;

  for (size_t j = q_size; j-- > 0;) {
    // D3: Compute trial quotient
    const DWord q_hat = computeTrialQuotient(u[j + n], u[j + n - 1], u[j + n - 2], v_high, v_mid);

    // D4: Multiply and subtract
    const SignedDWord borrow = multiplyAndSubtract(u, v, q_hat, j, n);

    // D5-D6: Store quotient digit, correct if needed
    q[j] = static_cast<Word>(q_hat);
    if (borrow < 0) {
      --q[j];
      addBack(u, v, j, n);
    }
  }

  // Build quotient
  BigInt quotient;
  quotient.digits_ = std::move(q);
  quotient.deleteZeroHighOrderDigit();
  quotient.positive_ = quotient_positive || !quotient;

  // D8: Unnormalize remainder
  BigInt remainder;
  remainder.digits_.resize(n);
  shiftWordsRight(u, remainder.digits_, n, shift);
  remainder.deleteZeroHighOrderDigit();
  remainder.positive_ = positive_ || !remainder;

  return {std::move(quotient), std::move(remainder)};
}

BigInt BigInt::operator/(const BigInt& divisor) const {
  return divmod(divisor).first;
}

BigInt& BigInt::operator/=(const BigInt& divisor) {
  *this = divmod(divisor).first;
  return *this;
}

BigInt BigInt::operator%(const BigInt& divisor) const {
  return divmod(divisor).second;
}

uint32_t BigInt::operator%(const uint32_t divisor) const {
  if (divisor == 0) {
    throw std::domain_error("Division by zero");
  }
  // Direct zero check to avoid template recursion
  if (digits_.empty() || (digits_.size() == 1 && digits_[0] == 0)) {
    return 0;
  }

  // Horner's method: ((d[n-1] * B + d[n-2]) * B + ...) mod divisor
  // where B = 2^32. We use the property that (a*B + b) mod m = ((a mod m)*B + b) mod m
  // Since remainder < divisor < 2^32 and B = 2^32, remainder*B fits in 64 bits.
  constexpr uint64_t base = kDigitBase;
  uint64_t remainder = 0;

  for (const auto& digit : std::views::reverse(digits_)) {
    // remainder is already < divisor, so remainder * base fits in 64 bits
    remainder = (remainder * base + digit) % divisor;
  }

  return static_cast<uint32_t>(remainder);
}

BigInt& BigInt::operator%=(const BigInt& divisor) {
  *this = divmod(divisor).second;
  return *this;
}

// ============================================================================
// Mathematical Functions
// ============================================================================

BigInt pow(const BigInt& base, const BigInt& exponent) {
  if (!exponent.positive_) {
    return BigInt(0);
  }

  BigInt result(1);
  result.digits_.reserve(base.digits_.size() * static_cast<size_t>(exponent));

  for (size_t bit_idx = exponent.bitLength() - 1; bit_idx > 0; --bit_idx) {
    if (exponent.digits_[bit_idx >> kDigitShift] & (1U << (bit_idx & kBitIndexMask))) {
      result *= base;
    }
    result *= result;
  }

  if (exponent % 2 != 0) {
    result *= base;
  }

  if (!base.positive_) {
    result.positive_ = (exponent % 2 == 0);
  }

  return result;
}

size_t log2(const BigInt& value) noexcept {
  return value.bitLength() - 1;
}

BigInt powmod(const BigInt& base, const BigInt& exponent, const BigInt& divisor) {
  if (divisor == 0) {
    throw std::domain_error("Modular exponentiation with zero divisor");
  }

  // Handle edge cases
  if (exponent == 0) {
    return BigInt(1);
  }
  if (base == 0) {
    return BigInt(0);
  }

  // Algorithm selection criteria:
  // - Montgomery: odd modulus, >= 8 words (256 bits), exponent >= 16 bits
  // - Barrett: any modulus, >= 4 words (128 bits), exponent >= 64 bits
  // - Standard: fallback for small numbers
  const size_t mod_words = divisor.digits_.size();
  const size_t exp_bits = exponent.bitLength();
  const bool is_odd_modulus = divisor.digits_[0] & 1;

  const bool use_montgomery = is_odd_modulus && mod_words >= 8 && exp_bits >= 16;
  const bool use_barrett = !use_montgomery && mod_words >= 4 && exp_bits >= 64;

  if (use_montgomery) {
    // Montgomery CIOS (for large odd moduli)
    internal::MontgomeryContext mont(divisor.digits_);
    const size_t k = mont.wordCount();

    // Convert base to Montgomery form
    internal::MontgomeryContext::WordVec base_vec(k, 0);
    BigInt base_mod = base % divisor;
    std::ranges::copy(base_mod.digits_, base_vec.begin());
    auto base_mont = mont.toMontgomery(base_vec);

    // Initialize result to 1 in Montgomery form
    internal::MontgomeryContext::WordVec one_vec(k, 0);
    one_vec[0] = 1;
    auto result_mont = mont.toMontgomery(one_vec);

    // Square-and-multiply in Montgomery form (left-to-right binary method)
    internal::MontgomeryContext::WordVec temp(k);

    for (size_t i = 0; i < exp_bits; ++i) {
      const size_t bit_idx = exp_bits - 1 - i;

      // Square
      mont.square(result_mont, temp);
      std::swap(result_mont, temp);

      // Multiply if bit is set
      if (exponent.testBit(bit_idx)) {
        mont.multiply(result_mont, base_mont, temp);
        std::swap(result_mont, temp);
      }
    }

    // Convert back from Montgomery form
    auto result_vec = mont.fromMontgomery(result_mont);

    // Build result BigInt
    BigInt result;
    result.digits_ = std::move(result_vec);
    result.deleteZeroHighOrderDigit();
    return result;
  }

  if (use_barrett) {
    // Barrett reduction (for even moduli or medium-sized odd moduli)
    internal::BarrettContext barrett(divisor);
    BigInt result(1);
    BigInt b = barrett.reduce(base);

    for (size_t i = 0; i < exp_bits; ++i) {
      const size_t bit_idx = exp_bits - 1 - i;
      result = barrett.mulmod(result, result);
      if (exponent.testBit(bit_idx)) {
        result = barrett.mulmod(result, b);
      }
    }
    return result;
  }

  // Standard square-and-multiply (for small numbers)
  BigInt result(1);
  BigInt b = base % divisor;

  for (size_t i = 0; i < exp_bits; ++i) {
    const size_t bit_idx = exp_bits - 1 - i;
    result = (result * result) % divisor;
    if (exponent.testBit(bit_idx)) {
      result = (result * b) % divisor;
    }
  }
  return result;
}

BigInt inversemod(BigInt value, const BigInt& modulus) {
  if (modulus == 0) {
    throw std::domain_error("Modular inverse with zero divisor");
  }
  if (!isCoprime(value, modulus)) {
    throw std::domain_error("Modular inverse does not exist (numbers not coprime)");
  }

  BigInt mod_copy(modulus);
  BigInt quotient;
  BigInt x0(0);
  BigInt x1(1);
  BigInt temp;

  while (value > 1) {
    quotient = value / mod_copy;
    temp = mod_copy;
    mod_copy = value % mod_copy;
    value = temp;
    temp = x0;
    x0 = x1 - (quotient * x0);
    x1 = temp;
  }

  if (!x1.positive_) {
    x1 += modulus;
  }

  return x1;
}

bool congruencemod(const BigInt& a, const BigInt& b, const BigInt& modulus) {
  if (modulus == 0) {
    throw std::domain_error("Congruence check with zero divisor");
  }

  BigInt rem_a = a % modulus;
  BigInt rem_b = b % modulus;

  while (!rem_a.positive_) {
    rem_a += modulus;
  }
  while (rem_a > modulus) {
    rem_a -= modulus;
  }
  while (!rem_b.positive_) {
    rem_b += modulus;
  }
  while (rem_b > modulus) {
    rem_b -= modulus;
  }

  return rem_a == rem_b;
}

bool isCoprime(const BigInt& a, const BigInt& b) {
  return gcd(a, b) == 1;
}

int8_t symbolJacobi(BigInt a, BigInt n) {
  if (!isCoprime(a, n)) {
    return 0;
  }

  int8_t result = 1;
  BigInt temp;

  if (!a.positive_) {
    a.positive_ = true;
    if (n % 4 == 3) {
      result = -result;
    }
  }

  while (a) {
    // Factor out powers of 2 using trailingZeros() instead of loop
    const size_t twos_count = a.trailingZeros();
    if (twos_count > 0) {
      a >>= twos_count;
    }

    if (twos_count % 2) {
      uint32_t n_mod_8 = n % 8;
      if (n_mod_8 == 3 || n_mod_8 == 5) {
        result = -result;
      }
    }

    if (a % 4 == 3 && n % 4 == 3) {
      result = -result;
    }

    temp = a;
    a = n % temp;
    n = temp;
  }

  return result;
}

// ============================================================================
// Bitwise Operators
// ============================================================================

BigInt BigInt::operator~() const {
  return -*this - BigInt(1);
}

BigInt BigInt::operator&(const BigInt& rhs) const {
  if (positive_ && rhs.positive_) {
    BigInt result;
    result.digits_.reserve(std::min(digits_.size(), rhs.digits_.size()));

    auto left_it = digits_.cbegin();
    auto right_it = rhs.digits_.cbegin();

    while (left_it != digits_.cend() && right_it != rhs.digits_.cend()) {
      result.digits_.emplace_back(*left_it & *right_it);
      ++left_it;
      ++right_it;
    }

    result.positive_ = true;
    return result;
  } else if (!positive_ && !rhs.positive_) {
    return -(~*this | ~rhs) - BigInt(1);
  } else if (positive_ && !rhs.positive_) {
    return (*this | ~rhs) + rhs + BigInt(1);
  } else {  // !positive_ && rhs.positive_
    return (~*this | rhs) + *this + BigInt(1);
  }
}

BigInt& BigInt::operator&=(const BigInt& rhs) {
  *this = *this & rhs;
  return *this;
}

BigInt BigInt::operator|(const BigInt& rhs) const {
  if (positive_ && rhs.positive_) {
    BigInt result;
    const bool this_larger = (digits_.size() >= rhs.digits_.size());

    result.digits_.reserve(this_larger ? digits_.size() + 1 : rhs.digits_.size() + 1);

    auto larger_it = this_larger ? digits_.cbegin() : rhs.digits_.cbegin();
    auto smaller_it = this_larger ? rhs.digits_.cbegin() : digits_.cbegin();
    auto larger_end = this_larger ? digits_.cend() : rhs.digits_.cend();
    auto smaller_end = this_larger ? rhs.digits_.cend() : digits_.cend();

    while (smaller_it != smaller_end) {
      result.digits_.emplace_back(*larger_it | *smaller_it);
      ++larger_it;
      ++smaller_it;
    }

    while (larger_it != larger_end) {
      result.digits_.emplace_back(*larger_it);
      ++larger_it;
    }

    result.positive_ = true;
    return result;
  } else if (!positive_ && !rhs.positive_) {
    return -(~*this & ~rhs) - BigInt(1);
  } else if (positive_ && !rhs.positive_) {
    return (*this & ~rhs) + rhs;
  } else {  // !positive_ && rhs.positive_
    return (~*this & rhs) + *this;
  }
}

BigInt& BigInt::operator|=(const BigInt& rhs) {
  *this = *this | rhs;
  return *this;
}

BigInt BigInt::operator^(const BigInt& rhs) const {
  if (positive_ && rhs.positive_) {
    BigInt result;
    const bool this_larger = (digits_.size() >= rhs.digits_.size());

    result.digits_.reserve(this_larger ? digits_.size() + 1 : rhs.digits_.size() + 1);

    auto larger_it = this_larger ? digits_.cbegin() : rhs.digits_.cbegin();
    auto smaller_it = this_larger ? rhs.digits_.cbegin() : digits_.cbegin();
    auto larger_end = this_larger ? digits_.cend() : rhs.digits_.cend();
    auto smaller_end = this_larger ? rhs.digits_.cend() : digits_.cend();

    while (smaller_it != smaller_end) {
      result.digits_.emplace_back(*larger_it ^ *smaller_it);
      ++larger_it;
      ++smaller_it;
    }

    while (larger_it != larger_end) {
      result.digits_.emplace_back(*larger_it);
      ++larger_it;
    }

    result.positive_ = true;
    result.deleteZeroHighOrderDigit();
    return result;
  } else {
    return (*this | rhs) & (~*this | ~rhs);
  }
}

BigInt& BigInt::operator^=(const BigInt& rhs) {
  *this = *this ^ rhs;
  return *this;
}

// ============================================================================
// Shift Operators
// ============================================================================

BigInt BigInt::operator<<(size_t shift) const {
  if (!shift || digits_.empty() || !(*this)) {
    return *this;
  }

  const size_t digit_shift = shift / kBitsPerDigit;
  const size_t bit_shift = shift & kBitIndexMask;

  BigInt result;
  result.positive_ = positive_;
  result.digits_.reserve(digits_.size() + digit_shift + 1);
  result.digits_.resize(digit_shift, 0);

  if (bit_shift == 0) {
    result.digits_.insert(result.digits_.end(), digits_.begin(), digits_.end());
  } else {
    uint32_t carry = 0;
    for (const auto& digit : digits_) {
      uint32_t shifted_digit = (digit << bit_shift) | carry;
      carry = digit >> (kBitsPerDigit - bit_shift);
      result.digits_.emplace_back(shifted_digit);
    }
    if (carry) {
      result.digits_.emplace_back(carry);
    }
  }

  return result;
}

BigInt& BigInt::operator<<=(size_t shift) {
  *this = *this << shift;
  return *this;
}

BigInt BigInt::operator>>(size_t shift) const {
  if (!shift || digits_.empty() || !(*this)) {
    return *this;
  }

  const size_t digit_shift = shift / kBitsPerDigit;
  const size_t bit_shift = shift & kBitIndexMask;

  if (digit_shift >= digits_.size()) {
    return BigInt(0);
  }

  BigInt result;
  result.positive_ = positive_;

  const size_t remaining_digits = digits_.size() - digit_shift;
  result.digits_.reserve(remaining_digits);

  if (bit_shift == 0) {
    result.digits_.insert(result.digits_.end(), digits_.begin() + static_cast<long>(digit_shift),
                          digits_.end());
  } else {
    const uint32_t mask = static_cast<uint32_t>((1ULL << bit_shift) - 1);
    uint32_t carry = 0;

    for (auto it = digits_.crbegin(); it != digits_.crend() - static_cast<long>(digit_shift);
         ++it) {
      uint32_t shifted_digit = (*it >> bit_shift) | carry;
      carry = (*it & mask) << (kBitsPerDigit - bit_shift);
      result.digits_.emplace_back(shifted_digit);
    }

    std::reverse(result.digits_.begin(), result.digits_.end());
  }

  result.deleteZeroHighOrderDigit();
  return result;
}

BigInt& BigInt::operator>>=(size_t shift) {
  *this = *this >> shift;
  return *this;
}

BigInt BigInt::leftCircularShift(size_t shift) const {
  const BigInt mask(--(BigInt(1) << bitLength()));
  return (((*this << shift) | (*this >> (bitLength() - shift))) & mask);
}

BigInt BigInt::rightCircularShift(size_t shift) const {
  const BigInt mask(--(BigInt(1) << bitLength()));
  return (((*this >> shift) | (*this << (bitLength() - shift))) & mask);
}

// ============================================================================
// Comparison Operators
// ============================================================================

bool BigInt::operator!() const noexcept {
  for (const auto& digit : digits_) {
    if (digit != 0) {
      return false;
    }
  }
  return true;
}

bool BigInt::operator&&(const BigInt& rhs) const noexcept {
  return static_cast<bool>(*this) && static_cast<bool>(rhs);
}

bool BigInt::operator||(const BigInt& rhs) const noexcept {
  return static_cast<bool>(*this) || static_cast<bool>(rhs);
}

std::strong_ordering BigInt::operator<=>(const BigInt& rhs) const noexcept {
  // Direct zero check to avoid recursion (don't use == 0 here)
  const bool this_zero = digits_.empty() || (digits_.size() == 1 && digits_[0] == 0);
  const bool rhs_zero = rhs.digits_.empty() || (rhs.digits_.size() == 1 && rhs.digits_[0] == 0);

  if (this_zero && rhs_zero) {
    return std::strong_ordering::equal;
  }

  if (positive_ && !rhs.positive_) {
    return std::strong_ordering::greater;
  }
  if (!positive_ && rhs.positive_) {
    return std::strong_ordering::less;
  }

  const int mag_cmp = compareMagnitude(rhs);

  if (positive_) {
    if (mag_cmp > 0)
      return std::strong_ordering::greater;
    if (mag_cmp < 0)
      return std::strong_ordering::less;
  } else {
    if (mag_cmp > 0)
      return std::strong_ordering::less;
    if (mag_cmp < 0)
      return std::strong_ordering::greater;
  }

  return std::strong_ordering::equal;
}

bool BigInt::operator==(const BigInt& rhs) const noexcept {
  return (*this <=> rhs) == std::strong_ordering::equal;
}

std::strong_ordering BigInt::compareToSigned(const int64_t rhs) const noexcept {
  // Direct zero check to avoid recursion (don't use == 0 here)
  const bool this_zero = digits_.empty() || (digits_.size() == 1 && digits_[0] == 0);

  // Fast path for zero comparison
  if (rhs == 0) {
    if (this_zero) {
      return std::strong_ordering::equal;
    }
    return positive_ ? std::strong_ordering::greater : std::strong_ordering::less;
  }

  // Sign comparison
  const bool rhs_positive = rhs > 0;
  if (positive_ != rhs_positive) {
    return positive_ ? std::strong_ordering::greater : std::strong_ordering::less;
  }

  // Both same sign - compare magnitudes
  const uint64_t rhs_abs = rhs > 0 ? static_cast<uint64_t>(rhs) : static_cast<uint64_t>(-rhs);

  // If we have more than 2 digits, we're definitely larger in magnitude
  if (digits_.size() > 2) {
    return positive_ ? std::strong_ordering::greater : std::strong_ordering::less;
  }

  // Get our magnitude as uint64_t
  uint64_t this_abs = 0;
  if (!digits_.empty()) {
    this_abs = digits_[0];
    if (digits_.size() > 1) {
      this_abs |= static_cast<uint64_t>(digits_[1]) << kBitsPerDigit;
    }
  }

  // Compare magnitudes, accounting for sign
  if (this_abs == rhs_abs) {
    return std::strong_ordering::equal;
  }
  if (this_abs > rhs_abs) {
    return positive_ ? std::strong_ordering::greater : std::strong_ordering::less;
  }
  return positive_ ? std::strong_ordering::less : std::strong_ordering::greater;
}

std::strong_ordering BigInt::compareToUnsigned(const uint64_t rhs) const noexcept {
  // Direct zero check to avoid recursion (don't use == 0 here)
  const bool this_zero = digits_.empty() || (digits_.size() == 1 && digits_[0] == 0);

  // Negative BigInt is always less than unsigned
  if (!positive_ && !this_zero) {
    return std::strong_ordering::less;
  }

  // Fast path for zero comparison
  if (rhs == 0) {
    return this_zero ? std::strong_ordering::equal : std::strong_ordering::greater;
  }

  // If we have more than 2 digits, we're definitely larger
  if (digits_.size() > 2) {
    return std::strong_ordering::greater;
  }

  // Get our magnitude as uint64_t
  uint64_t this_val = 0;
  if (!digits_.empty()) {
    this_val = digits_[0];
    if (digits_.size() > 1) {
      this_val |= static_cast<uint64_t>(digits_[1]) << kBitsPerDigit;
    }
  }

  if (this_val == rhs) {
    return std::strong_ordering::equal;
  }
  return this_val > rhs ? std::strong_ordering::greater : std::strong_ordering::less;
}

int BigInt::compareMagnitude(const BigInt& other) const noexcept {
  if (digits_.size() > other.digits_.size()) {
    return 1;
  } else if (digits_.size() < other.digits_.size()) {
    return -1;
  }

  // Compare from most significant digit (end of vector) to least significant
  for (size_t i = digits_.size(); i > 0; --i) {
    const auto& left_digit = digits_[i - 1];
    const auto& right_digit = other.digits_[i - 1];
    if (left_digit > right_digit) {
      return 1;
    }
    if (left_digit < right_digit) {
      return -1;
    }
  }

  return 0;
}

BigInt abs(const BigInt& value) {
  BigInt result(value);
  result.positive_ = true;
  return result;
}

BigInt sqrt(const BigInt& value) {
  if (value == 0) {
    return BigInt();
  }
  if (!value.positive_) {
    throw std::domain_error("Square root of negative number is undefined");
  }
  if (value == 1) {
    return BigInt(1);
  }

  size_t bit_len = value.bitLength();
  BigInt x = BigInt(1) << ((bit_len + 1) / 2);

  while (true) {
    BigInt x_new = (x + value / x) >> 1;
    if (x_new >= x) {
      break;
    }
    x = x_new;
  }

  return x;
}

/// @brief Computes GCD using Euclidean algorithm with Knuth Algorithm D division.
/// @details This implementation uses the classic Euclidean algorithm, which is
/// well-suited for BigInt because our optimized Knuth Algorithm D division
/// reduces operand size exponentially per iteration.
///
/// Alternative algorithms considered:
/// - Binary GCD (Stein's): Benchmarked ~6x slower due to more iterations
/// - Lehmer's algorithm: Complex, error-prone, marginal gains for readable code
BigInt gcd(BigInt a, BigInt b) {
  // Handle zero cases
  if (a == 0) {
    b.positive_ = true;
    return b;
  }
  if (b == 0) {
    a.positive_ = true;
    return a;
  }

  // Work with absolute values
  a.positive_ = true;
  b.positive_ = true;

  // Euclidean algorithm: gcd(a, b) = gcd(b, a mod b)
  while (b != 0) {
    BigInt r = a % b;
    a = std::move(b);
    b = std::move(r);
  }

  return a;
}

BigInt lcm(const BigInt& a, const BigInt& b) {
  return (a * b) / gcd(a, b);
}

const BigInt& max(const BigInt& a, const BigInt& b) noexcept {
  return a > b ? a : b;
}

const BigInt& min(const BigInt& a, const BigInt& b) noexcept {
  return a < b ? a : b;
}

void swap(BigInt& lhs, BigInt& rhs) noexcept {
  std::swap(lhs.positive_, rhs.positive_);
  std::swap(lhs.digits_, rhs.digits_);
}

// ============================================================================
// Number Theory
// ============================================================================

namespace {

/// @brief Miller-Rabin witness test helper.
bool millerRabinWitness(const BigInt& witness, const BigInt& d, size_t s, const BigInt& n) {
  const BigInt n_minus_1 = n - BigInt(1);
  BigInt x = powmod(witness, d, n);

  if (x == 1 || x == n_minus_1) {
    return true;
  }

  for (size_t r = 1; r < s; ++r) {
    x = powmod(x, BigInt(2), n);
    if (x == n_minus_1) {
      return true;
    }
    if (x == 1) {
      return false;
    }
  }

  return false;
}

}  // anonymous namespace

bool BigInt::isProbablePrime(size_t rounds) const {
  if (!positive_ || *this == 0) {
    return false;
  }
  if (*this == 2 || *this == 3) {
    return true;
  }
  if (*this % 2 == 0 || *this < 2) {
    return false;
  }

  // Quick divisibility check using small primes from anonymous namespace
  for (uint32_t p : kSmallPrimes) {
    if (*this == p) {
      return true;
    }
    if (*this % p == 0) {
      return false;
    }
  }

  // Write n-1 as 2^s * d where d is odd (using trailingZeros() instead of loop)
  BigInt d = *this - BigInt(1);
  const size_t s = d.trailingZeros();
  d >>= s;

  // Deterministic witnesses for small numbers
  if (*this < 2047) {
    return millerRabinWitness(BigInt(2), d, s, *this);
  }
  if (*this < 1373653) {
    return millerRabinWitness(BigInt(2), d, s, *this) && millerRabinWitness(BigInt(3), d, s, *this);
  }
  if (*this < 25326001) {
    return millerRabinWitness(BigInt(2), d, s, *this) &&
           millerRabinWitness(BigInt(3), d, s, *this) && millerRabinWitness(BigInt(5), d, s, *this);
  }
  if (*this < 3215031751ULL) {
    return millerRabinWitness(BigInt(2), d, s, *this) &&
           millerRabinWitness(BigInt(3), d, s, *this) &&
           millerRabinWitness(BigInt(5), d, s, *this) && millerRabinWitness(BigInt(7), d, s, *this);
  }

  // Probabilistic for larger numbers
  const BigInt n_minus_3 = *this - BigInt(3);
  for (size_t i = 0; i < rounds; ++i) {
    BigInt witness = BigInt::randomBelow(n_minus_3) + BigInt(2);
    if (!millerRabinWitness(witness, d, s, *this)) {
      return false;
    }
  }

  return true;
}

BigInt BigInt::randomBits(size_t num_bits) {
  if (num_bits == 0) {
    return BigInt();
  }

  static std::random_device rd;
  static std::mt19937_64 gen(rd());
  static std::uniform_int_distribution<uint32_t> dist32(0, UINT32_MAX);

  const size_t num_words = (num_bits + kBitIndexMask) / kBitsPerDigit;
  std::vector<uint32_t> digits(num_words);

  for (size_t i = 0; i < num_words; ++i) {
    digits[i] = dist32(gen);
  }

  size_t top_bits = num_bits & kBitIndexMask;
  if (top_bits == 0) {
    top_bits = kBitsPerDigit;
  }

  const uint32_t mask = (1U << top_bits) - 1;
  digits.back() &= mask;
  digits.back() |= (1U << (top_bits - 1));

  BigInt result;
  result.digits_ = std::move(digits);
  result.positive_ = true;
  return result;
}

BigInt BigInt::randomBelow(const BigInt& upper_bound) {
  if (upper_bound <= 1) {
    return BigInt();
  }

  size_t bits = upper_bound.bitLength();
  BigInt result;

  do {
    result = randomBits(bits);
  } while (result >= upper_bound);

  return result;
}

BigInt BigInt::randomPrime(size_t num_bits) {
  if (num_bits < 2) {
    return BigInt(2);
  }
  if (num_bits == 2) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<int> dist(0, 1);
    return BigInt(dist(gen) ? 2 : 3);
  }

  while (true) {
    BigInt candidate = randomBits(num_bits);

    if (candidate % 2 == 0) {
      candidate.digits_[0] |= 1;
    }

    bool divisible = false;
    for (uint32_t p : kSmallPrimes) {
      if (candidate == p) {
        return candidate;
      }
      if (candidate % p == 0) {
        divisible = true;
        break;
      }
    }

    if (divisible) {
      continue;
    }

    if (candidate.isProbablePrime()) {
      return candidate;
    }
  }
}

BigInt BigInt::nextPrime() const {
  if (*this <= 2) {
    return BigInt(2);
  }
  if (*this == 3) {
    return BigInt(3);
  }

  BigInt candidate = *this;
  if (candidate % 2 == 0) {
    candidate += BigInt(1);
  }

  const size_t bits = candidate.bitLength();
  const size_t max_iterations = std::max(static_cast<size_t>(1000000), bits * bits * 100);
  size_t iterations = 0;

  while (iterations < max_iterations) {
    ++iterations;

    bool divisible = false;
    for (uint32_t p : kSmallPrimes) {
      if (candidate == p) {
        return candidate;
      }
      if (candidate % p == 0) {
        divisible = true;
        break;
      }
    }

    if (!divisible && candidate.isProbablePrime()) {
      return candidate;
    }

    candidate += BigInt(2);
  }

  throw std::runtime_error("nextPrime: exceeded maximum iterations");
}

// ============================================================================
// I/O Operators
// ============================================================================

std::ostream& operator<<(std::ostream& out, const BigInt& value) {
  // Respect stream format flags like std::hex, std::uppercase, std::showbase
  const auto flags = out.flags();
  const auto base_field = flags & std::ios::basefield;

  // Determine base from stream flags
  uint8_t base = kBaseDecimal;
  if (base_field == std::ios::hex) {
    base = kBaseHexadecimal;
  } else if (base_field == std::ios::oct) {
    // Octal not supported, fall back to decimal
    base = kBaseDecimal;
  }

  // Handle zero specially
  if (value == 0) {
    if ((flags & std::ios::showbase) && base == kBaseHexadecimal) {
      out << (flags & std::ios::uppercase ? "0X0" : "0x0");
    } else {
      out << '0';
    }
    return out;
  }

  // Build the output string
  std::string str = value.toStdString(base);

  // Apply uppercase if requested (for hex)
  if ((flags & std::ios::uppercase) && base == kBaseHexadecimal) {
    for (char& c : str) {
      if (c >= 'a' && c <= 'f') {
        c = static_cast<char>(c - 'a' + 'A');
      }
    }
  }

  // Add base prefix if showbase is set
  if (flags & std::ios::showbase) {
    if (base == kBaseHexadecimal) {
      const bool is_negative = !str.empty() && str[0] == '-';
      const std::string prefix = (flags & std::ios::uppercase) ? "0X" : "0x";
      if (is_negative) {
        str.insert(1, prefix);  // Insert after minus sign
      } else {
        str.insert(0, prefix);
      }
    }
  }

  out << str;
  return out;
}

std::istream& operator>>(std::istream& in, BigInt& value) {
  std::string str;
  in >> str;
  value = BigInt(str);
  return in;
}

// ============================================================================
// Conversion Functions
// ============================================================================

std::string BigInt::toStdString(const uint8_t base) const {
  std::stringstream ss;

  if (digits_.empty() || !(*this)) {
    return "0";
  }

  if (!positive_) {
    ss << '-';
  }

  if (base == kBaseBinary) {
    for (const auto& digit : std::views::reverse(digits_)) {
      ss << std::bitset<sizeof(uint32_t) * 8>(digit);
    }
  } else if (base == kBaseHexadecimal) {
    for (const auto& digit : std::views::reverse(digits_)) {
      ss << std::hex << std::setw(8) << std::setfill('0') << digit;
    }
  } else {  // base == kBaseDecimal
    const BigInt decimal_repr = toBigIntDec();
    for (const auto& digit : std::views::reverse(decimal_repr.digits_)) {
      ss << std::dec << std::setw(9) << std::setfill('0') << digit;
    }
  }

  std::string result = ss.str();
  size_t start_pos = positive_ ? 0 : 1;
  size_t first_non_zero = result.find_first_not_of("-0");
  if (first_non_zero != std::string::npos && first_non_zero > start_pos) {
    result.erase(start_pos, first_non_zero - start_pos);
  }

  return result;
}

BigInt::operator std::vector<uint8_t>() const {
  if (digits_.empty()) {
    return {};  // Zero returns empty vector
  }

  std::vector<uint8_t> result;
  const size_t num_bytes = byteLength();
  result.reserve(num_bytes);

  const size_t partial_bytes = num_bytes % sizeof(uint32_t);

  // Handle partial bytes from the most significant digit
  if (partial_bytes > 0) {
    const uint32_t high_digit = digits_.back();
    for (size_t i = 0; i < partial_bytes; ++i) {
      result.emplace_back(static_cast<uint8_t>(high_digit >> ((partial_bytes - i - 1) * 8)));
    }
  }

  // Handle remaining full 32-bit digits (from second-to-last down to first)
  const size_t start_idx = partial_bytes > 0 ? digits_.size() - 1 : digits_.size();
  for (size_t idx = start_idx; idx > 0; --idx) {
    const uint32_t digit = digits_[idx - 1];
    for (uint8_t byte_idx = 0; byte_idx < sizeof(uint32_t); ++byte_idx) {
      result.emplace_back(static_cast<uint8_t>(digit >> ((sizeof(uint32_t) - byte_idx - 1) * 8)));
    }
  }

  return result;
}

uint64_t BigInt::toUint64() const noexcept {
  if (digits_.empty()) {
    return 0;
  }
  if (digits_.size() >= 2) {
    return (static_cast<uint64_t>(*std::next(digits_.cbegin())) << kBitsPerDigit) |
           static_cast<uint64_t>(*digits_.cbegin());
  }
  return digits_.front();
}

int64_t BigInt::toInt64() const noexcept {
  const auto magnitude = toUint64();
  return positive_ ? static_cast<int64_t>(magnitude) : -static_cast<int64_t>(magnitude);
}

BigInt::operator bool() const noexcept {
  for (const auto& digit : digits_) {
    if (digit != 0) {
      return true;
    }
  }
  return false;
}

size_t BigInt::bitLength() const noexcept {
  if (!(*this)) {
    return 1;
  }

  size_t len = (digits_.size() - 1) * sizeof(uint32_t) * 8;
  uint32_t high_digit = digits_.back();
  uint8_t high_bits = 0;

  while (high_digit) {
    high_digit >>= 1;
    ++high_bits;
  }

  return len + high_bits;
}

size_t BigInt::byteLength() const noexcept {
  if (!(*this)) {
    return 1;
  }

  size_t len = (digits_.size() - 1) * sizeof(uint32_t);
  uint32_t high_digit = digits_.back();
  uint8_t high_bytes = 0;

  while (high_digit) {
    high_digit >>= 8;
    ++high_bytes;
  }

  return len + high_bytes;
}

bool BigInt::testBit(size_t n) const noexcept {
  const size_t word_idx = n >> kDigitShift;
  if (word_idx >= digits_.size()) {
    return false;
  }
  const size_t bit_idx = n & kBitIndexMask;
  return (digits_[word_idx] >> bit_idx) & 1U;
}

size_t BigInt::trailingZeros() const noexcept {
  // Zero has no trailing zeros by convention
  if (digits_.size() == 1 && digits_[0] == 0) {
    return 0;
  }

  size_t count = 0;
  for (size_t i = 0; i < digits_.size(); ++i) {
    if (digits_[i] == 0) {
      count += kBitsPerDigit;
    } else {
      // Use compiler intrinsic for trailing zeros in a word
#if defined(_MSC_VER)
      unsigned long idx;
      _BitScanForward(&idx, digits_[i]);
      count += idx;
#else
      count += static_cast<size_t>(__builtin_ctz(digits_[i]));
#endif
      break;
    }
  }
  return count;
}

// ============================================================================
// Private Helper Methods
// ============================================================================

void BigInt::alignTo(BigInt& other) {
  if (digits_.size() > other.digits_.size()) {
    other.digits_.reserve(digits_.size());
    other.digits_.resize(digits_.size(), 0);
  } else if (other.digits_.size() > digits_.size()) {
    digits_.reserve(other.digits_.size());
    digits_.resize(other.digits_.size(), 0);
  }
}

void BigInt::deleteZeroHighOrderDigit() {
  while (digits_.size() > 1 && !digits_.back()) {
    digits_.pop_back();
  }
  if (digits_.empty()) {
    digits_.push_back(0);
  }
}

BigInt BigInt::shiftDigitsToHigh(size_t shift) const {
  BigInt result = *this;
  result.digits_.insert(result.digits_.begin(), shift, 0);
  return result;
}

BigInt BigInt::shiftDigitsToLow(size_t shift) const {
  BigInt result = *this;

  if (result.digits_.size() > shift) {
    result.digits_.erase(result.digits_.begin(), result.digits_.begin() + static_cast<long>(shift));
  } else {
    result.digits_.shrink_to_fit();
    result.digits_.reserve(1);
    result.digits_.emplace_back(0);
    result.positive_ = true;
  }

  return result;
}

/// @brief Schoolbook multiplication: O(n*m) row-oriented algorithm.
/// @details Row-by-row approach with running carry. Each result digit may be
/// updated multiple times, but cache locality is better than column-oriented
/// Comba for our use case. Benchmarked: row-oriented is ~15% faster than Comba.
BigInt BigInt::multiplySchoolbook(const BigInt& other) const {
  if (digits_.empty() || other.digits_.empty() || *this == 0 || other == 0) {
    return BigInt();
  }

  const size_t m = digits_.size();
  const size_t n = other.digits_.size();

  BigInt result;
  result.digits_.resize(m + n, 0);

  for (size_t i = 0; i < m; ++i) {
    uint64_t carry = 0;
    for (size_t j = 0; j < n; ++j) {
      uint64_t product =
          static_cast<uint64_t>(digits_[i]) * static_cast<uint64_t>(other.digits_[j]) +
          static_cast<uint64_t>(result.digits_[i + j]) + carry;
      result.digits_[i + j] = static_cast<uint32_t>(product);
      carry = product >> kBitsPerDigit;
    }
    result.digits_[i + n] = static_cast<uint32_t>(carry);
  }

  result.deleteZeroHighOrderDigit();
  result.positive_ = true;
  return result;
}

/// @brief Toom-Cook 3-way multiplication: O(n^1.465) complexity.
/// @details Splits each operand into 3 parts of size k = ceil(n/3):
///   p(x) = p0 + p1*x + p2*x^2  where x = B^k (B = 2^32)
///   q(x) = q0 + q1*x + q2*x^2
/// Evaluates at 5 points: 0, 1, -1, 2, ∞
/// Performs 5 recursive multiplications, then interpolates.
/// @note Intermediate values can be negative; signed arithmetic is used throughout.
BigInt BigInt::multiplyToom3(const BigInt& other) const {
  if (*this == 0 || other == 0) {
    return BigInt();
  }

  const size_t m = digits_.size();
  const size_t n = other.digits_.size();

  // Fall back to Karatsuba for smaller operands
  if (m < kToom3Threshold || n < kToom3Threshold) {
    return multiplyKaratsuba(other);
  }

  // Split into 3 parts of size k
  const size_t k = (std::max(m, n) + 2) / 3;

  // Helper: extract slice [start, start+len) from digits as BigInt
  auto sliceDigits = [](const std::vector<uint32_t>& digits, size_t start, size_t len) {
    BigInt result;
    if (start >= digits.size()) {
      return result;  // Zero
    }
    const auto begin = digits.begin() + static_cast<ptrdiff_t>(start);
    const auto end = begin + static_cast<ptrdiff_t>(std::min(len, digits.size() - start));
    result.digits_.assign(begin, end);
    result.deleteZeroHighOrderDigit();
    result.positive_ = true;
    return result;
  };

  // Split this = p0 + p1*B^k + p2*B^(2k)
  const BigInt p0 = sliceDigits(digits_, 0, k);
  const BigInt p1 = sliceDigits(digits_, k, k);
  const BigInt p2 = sliceDigits(digits_, 2 * k, k);

  // Split other = q0 + q1*B^k + q2*B^(2k)
  const BigInt q0 = sliceDigits(other.digits_, 0, k);
  const BigInt q1 = sliceDigits(other.digits_, k, k);
  const BigInt q2 = sliceDigits(other.digits_, 2 * k, k);

  // Helper: evaluate polynomial p0 + p1*x + p2*x^2 at point x
  auto evaluateAt = [](const BigInt& c0, const BigInt& c1, const BigInt& c2, int x) {
    switch (x) {
      case 0:
        return c0;
      case 1:
        return c0 + c1 + c2;
      case -1:
        return c0 - c1 + c2;  // Can be negative
      case 2:
        return c0 + (c1 << 1) + (c2 << 2);
      default:
        return c2;  // x = ∞ (leading coefficient)
    }
  };

  // Evaluate p(x) and q(x) at 5 points and multiply
  constexpr std::array kEvalPoints = {0, 1, -1, 2, 3};  // 3 represents ∞
  std::array<BigInt, 5> r{};
  for (size_t i = 0; i < kEvalPoints.size(); ++i) {
    const int x = kEvalPoints[i];
    r[i] = evaluateAt(p0, p1, p2, x) * evaluateAt(q0, q1, q2, x);
  }

  // Unpack evaluation results: r(0), r(1), r(-1), r(2), r(∞)
  const auto& [r_0, r_1, r_m1, r_2, r_inf] = r;

  // Interpolation to recover coefficients c0, c1, c2, c3, c4
  // r(x) = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4
  const BigInt c0 = r_0;
  const BigInt c4 = r_inf;

  // r_1 + r_m1 = 2*(c0 + c2 + c4), r_1 - r_m1 = 2*(c1 + c3)
  const BigInt sum_even = (r_1 + r_m1) / BigInt(2);  // c0 + c2 + c4
  const BigInt sum_odd = (r_1 - r_m1) / BigInt(2);   // c1 + c3
  const BigInt c2 = sum_even - c0 - c4;

  // Use r_2 to separate c1 and c3:
  // r_2 - c0 - 4*c2 - 16*c4 = 2*(c1 + 4*c3)
  const BigInt c1_4c3 = (r_2 - c0 - (c2 * BigInt(4)) - (c4 * BigInt(16))) / BigInt(2);
  const BigInt c3 = (c1_4c3 - sum_odd) / BigInt(3);
  const BigInt c1 = sum_odd - c3;

  // Combine: result = c0 + c1*B^k + c2*B^(2k) + c3*B^(3k) + c4*B^(4k)
  const std::array coeffs = {c0, c1, c2, c3, c4};
  BigInt result;
  for (size_t i = 0; i < coeffs.size(); ++i) {
    result = result + coeffs[i].shiftDigitsToHigh(i * k);
  }

  result.positive_ = true;
  return result;
}

BigInt BigInt::multiplyKaratsuba(const BigInt& other) const {
  if (*this == 0 || other == 0) {
    return BigInt();
  }

  const size_t m = digits_.size();
  const size_t n = other.digits_.size();

  if (m < kKaratsubaThreshold || n < kKaratsubaThreshold) {
    return multiplySchoolbook(other);
  }

  const size_t half = (std::max(m, n) + 1) / 2;

  // Split this = high1 * B^half + low1
  BigInt low1, high1;
  if (m <= half) {
    low1 = *this;
    low1.positive_ = true;
  } else {
    low1.digits_.assign(digits_.begin(), digits_.begin() + static_cast<long>(half));
    low1.positive_ = true;
    low1.deleteZeroHighOrderDigit();
    high1.digits_.assign(digits_.begin() + static_cast<long>(half), digits_.end());
    high1.positive_ = true;
    high1.deleteZeroHighOrderDigit();
  }

  // Split other = high2 * B^half + low2
  BigInt low2, high2;
  if (n <= half) {
    low2 = other;
    low2.positive_ = true;
  } else {
    low2.digits_.assign(other.digits_.begin(), other.digits_.begin() + static_cast<long>(half));
    low2.positive_ = true;
    low2.deleteZeroHighOrderDigit();
    high2.digits_.assign(other.digits_.begin() + static_cast<long>(half), other.digits_.end());
    high2.positive_ = true;
    high2.deleteZeroHighOrderDigit();
  }

  // Karatsuba: 3 multiplications instead of 4
  BigInt z0 = low1.multiplyKaratsuba(low2);
  BigInt z2 = high1.multiplyKaratsuba(high2);

  BigInt sum1 = abs(low1) + abs(high1);
  BigInt sum2 = abs(low2) + abs(high2);
  BigInt z1 = sum1.multiplyKaratsuba(sum2) - z0 - z2;

  BigInt result = z2.shiftDigitsToHigh(2 * half) + z1.shiftDigitsToHigh(half) + z0;
  result.positive_ = true;
  return result;
}

BigInt BigInt::toBigIntDec() const {
  const BigInt decimal_divisor(kDecimalBase);
  BigInt value = abs(*this);
  BigInt result;
  result.positive_ = positive_;
  result.digits_.reserve(digits_.size() + 1);

  while (value) {
    auto [quotient, remainder] = value.divmod(decimal_divisor);
    result.digits_.emplace_back(remainder.digits_.front());
    value = quotient;
  }

  return result;
}

}  // namespace bigint
