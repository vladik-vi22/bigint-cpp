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
 * - **Modular Exponentiation**: Montgomery CIOS for large odd moduli, standard square-and-multiply otherwise
 * - **Primality**: Miller-Rabin with deterministic witnesses for small numbers
 * - **GCD**: Euclidean algorithm with fast division
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
#include <bitset>
#include <cmath>
#include <compare>
#include <iomanip>
#include <iterator>
#include <random>
#include <ranges>
#include <sstream>
#include <stdexcept>

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
// Montgomery Arithmetic Helper (used exclusively by powmod)
// ============================================================================

/**
 * @brief Montgomery multiplication context for modular exponentiation.
 *
 * @details Implements the CIOS (Coarsely Integrated Operand Scanning) algorithm
 * which combines multiplication and reduction in a single pass, avoiding
 * temporary BigInt allocations.
 *
 * Montgomery form: x̃ = x * R mod N, where R = 2^(32*k) and k = word count of N.
 * Multiplication in Montgomery form: montMul(ã, b̃) = a*b*R mod N
 *
 * @note This struct is used exclusively by the `powmod` function for efficient
 *       modular exponentiation with large odd moduli. Montgomery multiplication
 *       replaces expensive division with cheaper addition/subtraction operations.
 *
 * @see "Handbook of Applied Cryptography", Chapter 14
 * @see https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
 */
struct MontgomeryContext {
  /// 32-bit word type for internal representation.
  using Word = uint32_t;
  /// 64-bit double-word type for intermediate calculations.
  using DWord = uint64_t;
  /// Vector of words for multi-precision arithmetic.
  using WordVec = std::vector<Word>;

  /**
   * @brief Constructs Montgomery context for given modulus.
   * @param modulus Odd modulus N (must be odd for Montgomery to work).
   * @pre modulus must be odd.
   */
  explicit MontgomeryContext(const WordVec& modulus)
      : n_(modulus), k_(modulus.size()), scratch_(k_ + 2, 0) {
    n0_inv_ = computeNegInverse(n_[0]);
    r_squared_ = computeRSquared();
  }

  /// Number of words in the modulus.
  [[nodiscard]] size_t wordCount() const noexcept { return k_; }

  /**
   * @brief Converts a value to Montgomery form: x̃ = x * R mod N.
   * @param x Input value (must be < N, zero-padded to k words).
   * @return x in Montgomery form.
   */
  [[nodiscard]] WordVec toMontgomery(const WordVec& x) const {
    WordVec result(k_);
    multiply(x, r_squared_, result);
    return result;
  }

  /**
   * @brief Converts from Montgomery form back to normal: x = x̃ * R^(-1) mod N.
   * @param x_mont Value in Montgomery form.
   * @return Normal representation.
   */
  [[nodiscard]] WordVec fromMontgomery(const WordVec& x_mont) const {
    WordVec one(k_, 0);
    one[0] = 1;
    WordVec result(k_);
    multiply(x_mont, one, result);
    return result;
  }

  /**
   * @brief Montgomery multiplication: result = a * b * R^(-1) mod N.
   * @param a First operand in Montgomery form.
   * @param b Second operand in Montgomery form.
   * @param[out] result Product in Montgomery form.
   */
  void multiply(const WordVec& a, const WordVec& b, WordVec& result) const {
    auto& t = scratch_;
    std::ranges::fill(t, Word{0});

    for (size_t i = 0; i < k_; ++i) {
      // Step 1: t += a[i] * b
      DWord carry = 0;
      for (size_t j = 0; j < k_; ++j) {
        DWord product = static_cast<DWord>(a[i]) * b[j] + t[j] + carry;
        t[j] = static_cast<Word>(product);
        carry = product >> kBitsPerDigit;
      }
      DWord sum = static_cast<DWord>(t[k_]) + carry;
      t[k_] = static_cast<Word>(sum);
      t[k_ + 1] = static_cast<Word>(sum >> kBitsPerDigit);

      // Step 2: m = t[0] * (-N^(-1)) mod 2^32
      Word m = t[0] * n0_inv_;

      // Step 3: t = (t + m * N) >> 32  (divide by word base)
      DWord mn0 = static_cast<DWord>(m) * n_[0] + t[0];
      carry = mn0 >> kBitsPerDigit;
      for (size_t j = 1; j < k_; ++j) {
        DWord product = static_cast<DWord>(m) * n_[j] + t[j] + carry;
        t[j - 1] = static_cast<Word>(product);
        carry = product >> kBitsPerDigit;
      }
      sum = static_cast<DWord>(t[k_]) + carry;
      t[k_ - 1] = static_cast<Word>(sum);
      t[k_] = t[k_ + 1] + static_cast<Word>(sum >> kBitsPerDigit);
      t[k_ + 1] = 0;
    }

    // Final reduction: if t >= N, subtract N
    conditionalSubtract(t, result);
  }

  /**
   * @brief Montgomery squaring: result = a^2 * R^(-1) mod N.
   * @param a Operand in Montgomery form.
   * @param[out] result Square in Montgomery form.
   * @note Currently delegates to multiply(). Could be optimized with dedicated
   *       squaring algorithm that exploits symmetry (a[i]*a[j] = a[j]*a[i]).
   */
  void square(const WordVec& a, WordVec& result) const { multiply(a, a, result); }

  /**
   * @brief Computes -N^(-1) mod 2^32 using Newton's method.
   * @param n0 Least significant word of N (must be odd).
   * @return The value -N^(-1) mod 2^32.
   * @details Newton iteration: x_{i+1} = x_i * (2 - n0 * x_i)
   *          Converges quadratically, reaching 32-bit precision in 5 iterations.
   */
  [[nodiscard]] static Word computeNegInverse(Word n0) noexcept {
    // Newton's method converges to n0^(-1) mod 2^32 in 5 iterations
    constexpr int kNewtonIterations = 5;
    Word inv = 1;
    for (int i = 0; i < kNewtonIterations; ++i) {
      inv *= 2 - n0 * inv;
    }
    return static_cast<Word>(-static_cast<int32_t>(inv));
  }

  /**
   * @brief Computes R^2 mod N for Montgomery conversion.
   * @return R^2 mod N as a word vector.
   */
  [[nodiscard]] WordVec computeRSquared() const {
    // R = 2^(32*k), so R^2 = 2^(64*k)
    // We compute this by shifting 1 left by 64*k bits, then taking mod N
    WordVec r2(2 * k_ + 1, 0);
    r2[2 * k_] = 1;  // R^2 = 2^(64*k)

    // Reduce R^2 mod N using repeated subtraction/shifting
    // This is only done once during setup, so efficiency is less critical
    return reduceModN(r2);
  }

  /**
   * @brief Converts little-endian word vector to BigInt.
   * @param words Little-endian word vector.
   * @return BigInt representation (always positive).
   */
  [[nodiscard]] static BigInt wordVecToBigInt(const WordVec& words) {
    // Reverse to big-endian for span constructor
    WordVec big_endian(words.rbegin(), words.rend());

    // Remove leading zeros using erase-remove idiom
    auto first_nonzero = std::ranges::find_if(big_endian, [](Word w) { return w != 0; });
    if (first_nonzero != big_endian.begin()) {
      big_endian.erase(big_endian.begin(), first_nonzero);
    }
    if (big_endian.empty()) {
      big_endian.push_back(0);
    }

    return BigInt(std::span<const Word>(big_endian), true);
  }

  /**
   * @brief Converts BigInt to little-endian word vector.
   * @param value BigInt to convert.
   * @param target_size Desired size of output vector (padded with zeros).
   * @return Little-endian word vector.
   */
  [[nodiscard]] static WordVec bigIntToWordVec(const BigInt& value, size_t target_size) {
    // Get bytes in big-endian order
    std::vector<uint8_t> bytes = static_cast<std::vector<uint8_t>>(value);

    // Convert bytes to words (big-endian bytes -> little-endian words)
    constexpr int kBitsPerByte = 8;
    WordVec result(target_size, 0);
    size_t byte_idx = bytes.size();
    for (size_t word_idx = 0; word_idx < target_size && byte_idx > 0; ++word_idx) {
      Word word = 0;
      for (size_t shift = 0; shift < kBitsPerDigit && byte_idx > 0; shift += kBitsPerByte) {
        word |= static_cast<Word>(bytes[--byte_idx]) << shift;
      }
      result[word_idx] = word;
    }
    return result;
  }

  /**
   * @brief Reduces a large value modulo N.
   * @param x Value to reduce (may be larger than N).
   * @return x mod N.
   * @note Uses BigInt's modulo operation internally. Only called during setup.
   */
  [[nodiscard]] WordVec reduceModN(const WordVec& x) const {
    if (x.empty() || std::ranges::all_of(x, [](Word w) { return w == 0; })) {
      return WordVec(k_, 0);
    }

    BigInt big_x = wordVecToBigInt(x);
    BigInt big_n = wordVecToBigInt(n_);
    BigInt result = big_x % big_n;

    return bigIntToWordVec(result, k_);
  }

  /**
   * @brief Conditionally subtracts N from t if t >= N.
   * @param t Input value (k+1 words, may be >= N).
   * @param[out] result Output value (k words, guaranteed < N).
   */
  void conditionalSubtract(const WordVec& t, WordVec& result) const {
    // Check if t >= N (compare from high to low)
    bool need_subtract = (t[k_] != 0);
    if (!need_subtract) {
      for (size_t i = k_; i-- > 0;) {
        if (t[i] > n_[i]) {
          need_subtract = true;
          break;
        }
        if (t[i] < n_[i]) {
          break;
        }
      }
    }

    if (need_subtract) {
      DWord borrow = 0;
      for (size_t i = 0; i < k_; ++i) {
        DWord diff = static_cast<DWord>(t[i]) - n_[i] - borrow;
        result[i] = static_cast<Word>(diff);
        borrow = (diff >> 63) & 1;
      }
    } else {
      std::ranges::copy_n(t.begin(), k_, result.begin());
    }
  }

  WordVec n_;                    ///< Modulus N
  size_t k_;                     ///< Word count of N
  Word n0_inv_;                  ///< -N^(-1) mod 2^32
  WordVec r_squared_;            ///< R^2 mod N for conversion to Montgomery form
  mutable WordVec scratch_;      ///< Scratch buffer for CIOS algorithm
};

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
  digits_.emplace_back(static_cast<uint32_t>(value & UINT32_MAX));
  digits_.emplace_back(static_cast<uint32_t>(value >> 32));
  deleteZeroHighOrderDigit();
}

void BigInt::initFromSigned(const int64_t value) {
  positive_ = value >= 0;
  const auto abs_value = static_cast<uint64_t>(value >= 0 ? value : -value);
  digits_.reserve(2);
  digits_.emplace_back(static_cast<uint32_t>(abs_value & UINT32_MAX));
  digits_.emplace_back(static_cast<uint32_t>(abs_value >> 32));
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
      sum.digits_.emplace_back(static_cast<uint32_t>(temp_sum & UINT32_MAX));
      carry = static_cast<uint32_t>(temp_sum >> 32);
      ++larger_it;
      ++smaller_it;
    }

    while (larger_it != larger_end) {
      temp_sum = static_cast<uint64_t>(*larger_it) + static_cast<uint64_t>(carry);
      sum.digits_.emplace_back(static_cast<uint32_t>(temp_sum & UINT32_MAX));
      carry = static_cast<uint32_t>(temp_sum >> 32);
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
      digits_[i] = static_cast<uint32_t>(sum & UINT32_MAX);
      carry = static_cast<uint32_t>(sum >> 32);
    }
    // Propagate carry through remaining digits
    for (; carry && i < digits_.size(); ++i) {
      uint64_t sum = static_cast<uint64_t>(digits_[i]) + carry;
      digits_[i] = static_cast<uint32_t>(sum & UINT32_MAX);
      carry = static_cast<uint32_t>(sum >> 32);
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
    product.digits_.emplace_back(static_cast<uint32_t>(temp_product & UINT32_MAX));
    carry = static_cast<uint32_t>(temp_product >> 32);
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
    digit = static_cast<uint32_t>(product & UINT32_MAX);
    carry = static_cast<uint32_t>(product >> 32);
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

  if (max_size >= kKaratsubaThreshold) {
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

std::pair<BigInt, BigInt> BigInt::DivMod(const BigInt& divisor) const {
  if (divisor == 0) {
    throw std::domain_error("Division by zero");
  }

  // Handle signs at the end - work with absolute values
  const int cmp = compareMagnitude(divisor);
  if (cmp < 0) {
    // |dividend| < |divisor|, quotient = 0, remainder = dividend
    BigInt quotient;
    BigInt remainder = *this;
    remainder.positive_ = true;
    if (!positive_ && remainder != 0) {
      remainder.positive_ = false;
    }
    return std::make_pair(quotient, remainder);
  }
  if (cmp == 0) {
    // |dividend| == |divisor|, quotient = ±1, remainder = 0
    BigInt quotient(1);
    quotient.positive_ = (positive_ == divisor.positive_);
    return std::make_pair(quotient, BigInt());
  }

  const size_t m = digits_.size();
  const size_t n = divisor.digits_.size();

  // Special case: single-word divisor - use fast path
  if (n == 1) {
    const uint32_t d = divisor.digits_[0];
    BigInt quotient;
    quotient.digits_.resize(m);
    uint64_t carry = 0;
    for (size_t i = m; i > 0; --i) {
      uint64_t cur = carry * kDigitBase + digits_[i - 1];
      quotient.digits_[i - 1] = static_cast<uint32_t>(cur / d);
      carry = cur % d;
    }
    quotient.deleteZeroHighOrderDigit();
    quotient.positive_ = (positive_ == divisor.positive_);
    if (quotient.digits_.empty()) {
      quotient.digits_.emplace_back(0);
      quotient.positive_ = true;
    }

    BigInt remainder(static_cast<uint32_t>(carry));
    if (!positive_ && remainder != 0) {
      remainder.positive_ = false;
    }
    return std::make_pair(quotient, remainder);
  }

  // Knuth's Algorithm D (TAOCP Vol 2, Section 4.3.1)
  // Divides u[0..m-1] by v[0..n-1] where m >= n >= 2

  // D1: Normalize - shift so that v[n-1] >= base/2
  const uint32_t v_high = divisor.digits_[n - 1];
  size_t shift = 0;
  if (v_high != 0) {
    // Count leading zeros using bit manipulation
    uint32_t temp = v_high;
    while ((temp & 0x80000000) == 0) {
      ++shift;
      temp <<= 1;
    }
  }

  // Create normalized copies
  std::vector<uint32_t> u(m + 1, 0);  // dividend with extra word
  std::vector<uint32_t> v(n, 0);       // divisor

  // Shift divisor left by 'shift' bits
  if (shift == 0) {
    for (size_t i = 0; i < n; ++i) {
      v[i] = divisor.digits_[i];
    }
  } else {
    uint32_t carry = 0;
    for (size_t i = 0; i < n; ++i) {
      uint32_t new_val = (divisor.digits_[i] << shift) | carry;
      carry = divisor.digits_[i] >> (kBitsPerDigit - shift);
      v[i] = new_val;
    }
  }

  // Shift dividend left by 'shift' bits
  if (shift == 0) {
    for (size_t i = 0; i < m; ++i) {
      u[i] = digits_[i];
    }
    u[m] = 0;
  } else {
    uint32_t carry = 0;
    for (size_t i = 0; i < m; ++i) {
      uint32_t new_val = (digits_[i] << shift) | carry;
      carry = digits_[i] >> (kBitsPerDigit - shift);
      u[i] = new_val;
    }
    u[m] = carry;
  }

  // D2-D7: Main loop - compute quotient digits from high to low
  const size_t q_size = m - n + 1;
  std::vector<uint32_t> q(q_size, 0);

  const uint64_t v_n1 = v[n - 1];  // Most significant word of divisor
  const uint64_t v_n2 = (n >= 2) ? v[n - 2] : 0;

  for (size_t j = q_size; j > 0; --j) {
    const size_t idx = j - 1;  // Current quotient position

    // D3: Calculate trial quotient q_hat
    const uint64_t u_high = (static_cast<uint64_t>(u[idx + n]) << kBitsPerDigit) + u[idx + n - 1];
    uint64_t q_hat = u_high / v_n1;
    uint64_t r_hat = u_high % v_n1;

    // Refine q_hat using second digit
    while (q_hat >= kDigitBase ||
           q_hat * v_n2 > (r_hat << kBitsPerDigit) + u[idx + n - 2]) {
      --q_hat;
      r_hat += v_n1;
      if (r_hat >= kDigitBase) {
        break;
      }
    }

    // D4: Multiply and subtract: u[idx..idx+n] -= q_hat * v[0..n-1]
    int64_t borrow = 0;
    for (size_t i = 0; i < n; ++i) {
      uint64_t product = q_hat * v[i];
      int64_t diff = static_cast<int64_t>(u[idx + i]) - static_cast<int64_t>(product & 0xFFFFFFFF) + borrow;
      u[idx + i] = static_cast<uint32_t>(diff & 0xFFFFFFFF);
      borrow = (diff >> kBitsPerDigit) - static_cast<int64_t>(product >> kBitsPerDigit);
    }
    int64_t diff = static_cast<int64_t>(u[idx + n]) + borrow;
    u[idx + n] = static_cast<uint32_t>(diff & 0xFFFFFFFF);

    // D5: Test remainder - if negative, q_hat was too large
    q[idx] = static_cast<uint32_t>(q_hat);
    if (diff < 0) {
      // D6: Add back - this happens rarely (probability ~2/base)
      --q[idx];
      uint64_t carry = 0;
      for (size_t i = 0; i < n; ++i) {
        uint64_t sum = static_cast<uint64_t>(u[idx + i]) + v[i] + carry;
        u[idx + i] = static_cast<uint32_t>(sum & 0xFFFFFFFF);
        carry = sum >> kBitsPerDigit;
      }
      u[idx + n] += static_cast<uint32_t>(carry);
    }
  }

  // Build quotient
  BigInt quotient;
  quotient.digits_ = std::move(q);
  quotient.deleteZeroHighOrderDigit();
  quotient.positive_ = (positive_ == divisor.positive_);
  if (quotient.digits_.empty()) {
    quotient.digits_.emplace_back(0);
    quotient.positive_ = true;
  }

  // D8: Unnormalize remainder - shift right by 'shift' bits
  BigInt remainder;
  remainder.digits_.resize(n);
  if (shift == 0) {
    for (size_t i = 0; i < n; ++i) {
      remainder.digits_[i] = u[i];
    }
  } else {
    uint32_t carry = 0;
    for (size_t i = n; i > 0; --i) {
      uint32_t new_val = (u[i - 1] >> shift) | carry;
      carry = u[i - 1] << (kBitsPerDigit - shift);
      remainder.digits_[i - 1] = new_val;
    }
  }
  remainder.deleteZeroHighOrderDigit();
  remainder.positive_ = true;
  if (remainder.digits_.empty()) {
    remainder.digits_.emplace_back(0);
  }

  // Apply sign to remainder (C++ truncated division: remainder has same sign as dividend)
  if (!positive_ && remainder != 0) {
    remainder.positive_ = false;
  }

  return std::make_pair(quotient, remainder);
}

BigInt BigInt::operator/(const BigInt& divisor) const {
  return DivMod(divisor).first;
}

BigInt& BigInt::operator/=(const BigInt& divisor) {
  *this = DivMod(divisor).first;
  return *this;
}

BigInt BigInt::operator%(const BigInt& divisor) const {
  return DivMod(divisor).second;
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
  const uint64_t base = static_cast<uint64_t>(1) << 32;
  uint64_t remainder = 0;

  for (const auto& digit : std::views::reverse(digits_)) {
    // remainder is already < divisor, so remainder * base fits in 64 bits
    remainder = (remainder * base + digit) % divisor;
  }

  return static_cast<uint32_t>(remainder);
}

BigInt& BigInt::operator%=(const BigInt& divisor) {
  *this = DivMod(divisor).second;
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
    if (exponent.digits_[bit_idx >> 5] & (1 << (bit_idx & 31))) {
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

  // Montgomery multiplication requires odd modulus and benefits from larger numbers.
  // For small moduli or even moduli, use standard square-and-multiply with fast division.
  // Threshold: Montgomery is beneficial when modulus has >= 8 words (256 bits) AND
  // exponent has >= 16 bits (enough operations to amortize setup cost).
  // Note: RSA public exponent 65537 is 17 bits, so we use 16 as threshold.
  const bool use_montgomery = (divisor.digits_[0] & 1) && divisor.digits_.size() >= 8 &&
                              exponent.bitLength() >= 16;

  if (!use_montgomery) {
    BigInt result(1);
    BigInt b = base % divisor;
    const size_t bit_len = exponent.bitLength();

    for (size_t i = 0; i < bit_len; ++i) {
      size_t bit_idx = bit_len - 1 - i;
      result = (result * result) % divisor;
      if (exponent.digits_[bit_idx >> 5] & (1 << (bit_idx & 31))) {
        result = (result * b) % divisor;
      }
    }
    return result;
  }

  // Montgomery exponentiation using CIOS algorithm
  MontgomeryContext mont(divisor.digits_);
  const size_t k = mont.wordCount();

  // Convert base to Montgomery form
  MontgomeryContext::WordVec base_vec(k, 0);
  BigInt base_mod = base % divisor;
  std::ranges::copy(base_mod.digits_, base_vec.begin());
  auto base_mont = mont.toMontgomery(base_vec);

  // Initialize result to 1 in Montgomery form
  MontgomeryContext::WordVec one_vec(k, 0);
  one_vec[0] = 1;
  auto result_mont = mont.toMontgomery(one_vec);

  // Square-and-multiply in Montgomery form (left-to-right binary method)
  MontgomeryContext::WordVec temp(k);
  const size_t bit_len = exponent.bitLength();

  for (size_t i = 0; i < bit_len; ++i) {
    const size_t bit_idx = bit_len - 1 - i;
    const bool bit_set = exponent.digits_[bit_idx >> kDigitShift] & (1U << (bit_idx & kBitIndexMask));

    // Square
    mont.square(result_mont, temp);
    std::swap(result_mont, temp);

    // Multiply if bit is set
    if (bit_set) {
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
    size_t twos_count = 0;
    while (a % 2 == 0) {
      a >>= 1;
      ++twos_count;
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
      this_abs |= static_cast<uint64_t>(digits_[1]) << 32;
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
      this_val |= static_cast<uint64_t>(digits_[1]) << 32;
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

  // Euclidean algorithm using fast division (Knuth Algorithm D)
  // Now that division is fast, this is efficient
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

  // Write n-1 as 2^s * d where d is odd
  BigInt d = *this - BigInt(1);
  size_t s = 0;
  while (d % 2 == 0) {
    d >>= 1;
    ++s;
  }

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
    return (static_cast<uint64_t>(*std::next(digits_.cbegin())) << 32) |
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
      result.digits_[i + j] = static_cast<uint32_t>(product & UINT32_MAX);
      carry = product >> 32;
    }
    result.digits_[i + n] = static_cast<uint32_t>(carry);
  }

  result.deleteZeroHighOrderDigit();
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
    auto [quotient, remainder] = value.DivMod(decimal_divisor);
    result.digits_.emplace_back(remainder.digits_.front());
    value = quotient;
  }

  return result;
}

}  // namespace bigint


