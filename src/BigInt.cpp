#include <bigint/BigInt.hpp>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <compare>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>

namespace bigint {

namespace {

// Internal constants (not exposed in public API)
constexpr uint64_t kBasisCalcSys = 1ULL << 32;  // 2^32 for carry calculations
constexpr uint32_t kBasisCalcDec = 1000000000;  // 10^9 for decimal conversion
constexpr uint8_t kDecimalCellSize = 9;         // Digits per decimal cell
constexpr size_t kKaratsubaThreshold = 32;      // Threshold for Karatsuba multiplication

// Small primes for quick divisibility rejection in prime generation
constexpr auto kSmallPrimes = std::to_array<uint32_t>({
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113
});

/// Converts a decimal string to a binary string.
/// Internal helper function for string constructor.
std::string DecimalToBinaryString(std::string decimal_str) {
  if (decimal_str == "0") {
    return "0";
  }
  std::string binary_str;
  std::vector<uint32_t> digits;
  std::vector<uint32_t> zero_arr;
  uint32_t carry_next;
  uint32_t carry_current;

  while (decimal_str.length() % kDecimalCellSize) {
    decimal_str.insert(0, 1, '0');
  }
  size_t num_cells = decimal_str.length() / kDecimalCellSize;
  digits.reserve(num_cells);
  for (size_t i = 0; i < num_cells; ++i) {
    digits.emplace_back(static_cast<uint32_t>(
        std::stoul(decimal_str.substr(i * kDecimalCellSize, kDecimalCellSize),
                   nullptr, 10)));
  }
  zero_arr.resize(digits.size(), 0);
  while (digits != zero_arr) {
    carry_next = 0;
    for (auto it = digits.begin(); it != digits.end(); ++it) {
      carry_current = carry_next;
      carry_next = (*it & 1);
      *it = (*it + carry_current * kBasisCalcDec) >> 1;
    }
    binary_str.insert(binary_str.begin(), 1, carry_next ? '1' : '0');
  }
  return binary_str;
}

}  // anonymous namespace



BigInt::BigInt() : positive_(true), digits_() {
}

BigInt::BigInt(const BigInt& other)
    : positive_(other.positive_), digits_(other.digits_) {
}

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

  const uint8_t cell_size = base == kBaseHexadecimal
                                ? (sizeof(uint32_t) * 2)
                                : (sizeof(uint32_t) * 8);
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
      throw std::invalid_argument(
          "Invalid character '" + std::string(1, str[i]) +
          "' at position " + std::to_string(i) +
          " for base " + std::to_string(base));
    }
  }

  if (base == kBaseDecimal) {
    str = DecimalToBinaryString(str);
  }
  while (str.length() % cell_size) {
    str.insert(0, 1, '0');
  }
  size_t num_cells = str.length() / cell_size;
  digits_.reserve(num_cells);
  for (size_t i = 0; i < num_cells; ++i) {
    digits_.emplace(
        digits_.begin(),
        static_cast<uint32_t>(std::stoul(
            str.substr(i * cell_size, cell_size), nullptr,
            base == kBaseHexadecimal ? kBaseHexadecimal : kBaseBinary)));
  }
}

BigInt::BigInt(const std::vector<uint32_t>& vec, const bool is_positive_)
    : positive_(is_positive_), digits_(vec) {
  std::reverse(digits_.begin(), digits_.end());
  deleteZeroHighOrderDigit();
}

BigInt::BigInt(const std::vector<uint16_t>& vec, const bool is_positive_) {
  digits_.reserve(vec.size() & 1 ? (vec.size() >> 1) + 1 : vec.size() >> 1);
  auto it = vec.crbegin();
  for (size_t i = 0; i < (vec.size() >> 1); ++i) {
    digits_.emplace_back(static_cast<uint32_t>(*it) |
                         static_cast<uint32_t>(*(++it)) << 16);
    ++it;
  }
  if (vec.size() & 1) {
    digits_.emplace_back(static_cast<uint32_t>(*it));
  }
  deleteZeroHighOrderDigit();
  positive_ = is_positive_;
}

BigInt::BigInt(const std::vector<uint8_t>& vec, const bool is_positive_) {
  digits_.reserve(vec.size() & 3 ? (vec.size() >> 2) + 1 : vec.size() >> 2);
  auto it = vec.crbegin();
  for (size_t i = 0; i < (vec.size() >> 2); ++i) {
    digits_.emplace_back(static_cast<uint32_t>(*it) |
                         static_cast<uint32_t>(*(++it)) << 8 |
                         static_cast<uint32_t>(*(++it)) << 16 |
                         static_cast<uint32_t>(*(++it)) << 24);
    ++it;
  }
  if ((vec.size() & 3) == 3) {
    digits_.emplace_back(static_cast<uint32_t>(*it) |
                         static_cast<uint32_t>(*(++it)) << 8 |
                         static_cast<uint32_t>(*(++it)) << 16);
  } else if ((vec.size() & 3) == 2) {
    digits_.emplace_back(static_cast<uint32_t>(*it) |
                         static_cast<uint32_t>(*(++it)) << 8);
  } else if ((vec.size() & 3) == 1) {
    digits_.emplace_back(static_cast<uint32_t>(*it));
  }
  deleteZeroHighOrderDigit();
  positive_ = is_positive_;
}

BigInt::BigInt(const std::vector<bool>& vec, const bool is_positive_) {
  digits_.reserve(vec.size() & 31 ? (vec.size() >> 5) + 1 : vec.size() >> 5);
  uint32_t element;
  auto it = vec.crbegin();
  for (size_t i = 0; i < (vec.size() >> 5); ++i) {
    element = 0;
    for (uint8_t bit_idx = 0; bit_idx < 32; ++bit_idx) {
      element |= static_cast<uint32_t>(*it) << bit_idx;
      ++it;
    }
    digits_.emplace_back(element);
  }
  if (vec.size() & 31) {
    element = 0;
    for (uint8_t bit_idx = 0; bit_idx < (vec.size() & 31); ++bit_idx) {
      element |= static_cast<uint32_t>(*it) << bit_idx;
      ++it;
    }
    digits_.emplace_back(element);
  }
  deleteZeroHighOrderDigit();
  positive_ = is_positive_;
}

BigInt::BigInt(const uint64_t value, const bool is_positive_)
    : positive_(is_positive_) {
  digits_.reserve(2);
  digits_.emplace_back(static_cast<uint32_t>(value & UINT32_MAX));
  digits_.emplace_back(static_cast<uint32_t>(value >> 32));
}

BigInt::BigInt(const uint32_t value, const bool is_positive_)
    : positive_(is_positive_), digits_{value} {
}

BigInt::BigInt(const int64_t value)
    : positive_(value >= 0) {
  const auto abs_value = static_cast<uint64_t>(std::abs(value));
  digits_.reserve(2);
  digits_.emplace_back(static_cast<uint32_t>(abs_value & UINT32_MAX));
  digits_.emplace_back(static_cast<uint32_t>(abs_value >> 32));
}

BigInt::BigInt(const int32_t value)
    : positive_(value >= 0), digits_{static_cast<uint32_t>(std::abs(value))} {
}

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

BigInt BigInt::operator+() const {
  return *this;
}

BigInt BigInt::operator+(const BigInt& addend) const {
  if (positive_ && addend.positive_) {
    BigInt sum;
    uint32_t carry = 0;
    uint64_t temp_sum;
    const bool this_larger = (digits_.size() >= addend.digits_.size());

    sum.digits_.reserve(this_larger ? digits_.size() + 1
                                    : addend.digits_.size() + 1);

    auto larger_it = this_larger ? digits_.cbegin() : addend.digits_.cbegin();
    auto smaller_it = this_larger ? addend.digits_.cbegin() : digits_.cbegin();
    auto larger_end = this_larger ? digits_.cend() : addend.digits_.cend();
    auto smaller_end = this_larger ? addend.digits_.cend() : digits_.cend();

    while (smaller_it != smaller_end) {
      temp_sum = static_cast<uint64_t>(*larger_it) +
                 static_cast<uint64_t>(*smaller_it) +
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
  if (addend.isZero()) {
    return *this;
  }
  if (isZero()) {
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
      uint64_t sum = static_cast<uint64_t>(digits_[i]) +
                     static_cast<uint64_t>(addend.digits_[i]) +
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
          (i < smaller->digits_.size())
              ? static_cast<int64_t>(smaller->digits_[i])
              : 0;
      int64_t diff =
          static_cast<int64_t>(larger->digits_[i]) - smaller_digit - borrow;
      if (diff >= 0) {
        result[i] = static_cast<uint32_t>(diff);
        borrow = 0;
      } else {
        result[i] = static_cast<uint32_t>(diff + static_cast<int64_t>(kBasisCalcSys));
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
        if (isZero()) {
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
                            static_cast<int64_t>(*subtrahend_it) -
                            static_cast<int64_t>(borrow);
        if (temp_diff >= 0) {
          diff.digits_.emplace_back(static_cast<uint32_t>(temp_diff));
          borrow = 0;
        } else {
          diff.digits_.emplace_back(
              static_cast<uint32_t>(temp_diff + static_cast<int64_t>(kBasisCalcSys)));
          borrow = 1;
        }
        ++minuend_it;
        ++subtrahend_it;
      }

      while (minuend_it != digits_.cend()) {
        int64_t temp_diff =
            static_cast<int64_t>(*minuend_it) - static_cast<int64_t>(borrow);
        if (temp_diff >= 0) {
          diff.digits_.emplace_back(static_cast<uint32_t>(temp_diff));
          borrow = 0;
        } else {
          diff.digits_.emplace_back(
              static_cast<uint32_t>(temp_diff + static_cast<int64_t>(kBasisCalcSys)));
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
  if (subtrahend.isZero()) {
    return *this;
  }
  BigInt negated = subtrahend;
  negated.positive_ = !subtrahend.positive_;
  return *this += negated;
}

BigInt& BigInt::operator--() {
  if (positive_) {
    if (digits_.empty() || isZero()) {
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

  for (auto it = digits_.cbegin(); it != digits_.cend(); ++it) {
    uint64_t temp_product =
        static_cast<uint64_t>(*it) * static_cast<uint64_t>(multiplier) +
        static_cast<uint64_t>(carry);
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
  if (isZero() || multiplier.isZero()) {
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
  if (divisor.isZero()) {
    throw std::domain_error("Division by zero");
  }

  const BigInt abs_divisor = abs(divisor);
  const size_t divisor_bit_len = abs_divisor.bitLength();

  BigInt quotient;
  quotient.digits_.reserve(digits_.size());

  BigInt remainder = abs(*this);

  if (remainder < abs_divisor) {
    quotient.positive_ = (positive_ == divisor.positive_);
    remainder.positive_ = positive_;
    return std::make_pair(quotient, remainder);
  }

  BigInt shifted_divisor;
  shifted_divisor.digits_.reserve(remainder.digits_.size());

  while (remainder >= abs_divisor) {
    size_t bit_diff = remainder.bitLength() - divisor_bit_len;
    shifted_divisor = abs_divisor << bit_diff;

    if (remainder < shifted_divisor) {
      shifted_divisor >>= 1;
      --bit_diff;
    }

    remainder -= shifted_divisor;

    // Create power of 2 directly
    BigInt power_of_two;
    size_t word_idx = bit_diff / 32;
    size_t bit_idx = bit_diff % 32;
    power_of_two.digits_.resize(word_idx + 1, 0);
    power_of_two.digits_[word_idx] = static_cast<uint32_t>(1) << bit_idx;
    quotient += power_of_two;
  }

  quotient.positive_ = (positive_ == divisor.positive_);
  remainder.positive_ = positive_;

  if (divisor.positive_) {
    while (!remainder.positive_) {
      remainder += divisor;
    }
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

BigInt& BigInt::operator%=(const BigInt& divisor) {
  *this = DivMod(divisor).second;
  return *this;
}

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

  if (exponent.isOdd()) {
    result *= base;
  }

  if (!base.positive_) {
    result.positive_ = exponent.isEven();
  }

  return result;
}

size_t log2(const BigInt& value) noexcept {
  return value.bitLength() - 1;
}

// Barrett reduction (commented out - kept for reference)
// BigInt powmod(BigInt base, const BigInt& exponent, const BigInt& divisor) {
//   BigInt power(1);
//   power.digits_.reserve(divisor.digits_.size());
//   const BigInt mu = power.shiftDigitsToHigh(divisor.digits_.size() * 2) / divisor;
//   const uint32_t bit_len = exponent.bitLength();
//   for (size_t bit_idx = 0; bit_idx < bit_len; ++bit_idx) {
//     if (exponent.digits_[bit_idx >> 5] & (1 << (bit_idx & 31))) {
//       power = BarrettReduction(power * base, divisor, mu);
//     }
//     base = BarrettReduction(base * base, divisor, mu);
//   }
//   return power;
// }

BigInt powmod(const BigInt& base, const BigInt& exponent, const BigInt& divisor) {
  if (divisor.isZero()) {
    throw std::domain_error("Modular exponentiation with zero divisor");
  }

  BigInt result(1);
  result.digits_.reserve(divisor.digits_.size());
  const size_t bit_len = exponent.bitLength();

  for (size_t i = 0; i < bit_len; ++i) {
    size_t bit_idx = bit_len - 1 - i;
    result = (result * result) % divisor;
    if (exponent.digits_[bit_idx >> 5] & (1 << (bit_idx & 31))) {
      result = (result * base) % divisor;
    }
  }

  return result;
}

BigInt inversemod(BigInt value, const BigInt& modulus) {
  if (modulus.isZero()) {
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

  while (value > BigInt(1)) {
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
  if (modulus.isZero()) {
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
  return gcd(a, b) == BigInt(1);
}

int8_t symbolJacobi(BigInt a, BigInt n) {
  if (!isCoprime(a, n)) {
    return 0;
  }

  int8_t result = 1;
  BigInt temp;

  if (!a.positive_) {
    a.positive_ = true;
    if (n % BigInt(4) == BigInt(3)) {
      result = -result;
    }
  }

  while (a) {
    size_t twos_count = 0;
    while (a.isEven()) {
      a >>= 1;
      ++twos_count;
    }

    if (twos_count % 2) {
      BigInt n_mod_8 = n % BigInt(8);
      if (n_mod_8 == BigInt(3) || n_mod_8 == BigInt(5)) {
        result = -result;
      }
    }

    if (a % BigInt(4) == BigInt(3) && n % BigInt(4) == BigInt(3)) {
      result = -result;
    }

    temp = a;
    a = n % temp;
    n = temp;
  }

  return result;
}

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

    result.digits_.reserve(this_larger ? digits_.size() + 1
                                       : rhs.digits_.size() + 1);

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

    result.digits_.reserve(this_larger ? digits_.size() + 1
                                       : rhs.digits_.size() + 1);

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
    return result;
  } else {
    return (*this | rhs) & (~*this | ~rhs);
  }
}

BigInt& BigInt::operator^=(const BigInt& rhs) {
  *this = *this ^ rhs;
  return *this;
}

BigInt BigInt::operator<<(size_t shift) const {
  if (!shift || digits_.empty() || !(*this)) {
    return *this;
  }

  const size_t digit_shift = shift / 32;
  const size_t bit_shift = shift % 32;

  BigInt result;
  result.positive_ = positive_;
  result.digits_.reserve(digits_.size() + digit_shift + 1);
  result.digits_.resize(digit_shift, 0);

  if (bit_shift == 0) {
    result.digits_.insert(result.digits_.end(), digits_.begin(), digits_.end());
  } else {
    uint32_t carry = 0;
    for (auto it = digits_.cbegin(); it != digits_.cend(); ++it) {
      uint32_t shifted_digit = (*it << bit_shift) | carry;
      carry = *it >> (32 - bit_shift);
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

  const size_t digit_shift = shift / 32;
  const size_t bit_shift = shift % 32;

  if (digit_shift >= digits_.size()) {
    return BigInt(0);
  }

  BigInt result;
  result.positive_ = positive_;

  const size_t remaining_digits = digits_.size() - digit_shift;
  result.digits_.reserve(remaining_digits);

  if (bit_shift == 0) {
    result.digits_.insert(result.digits_.end(),
                          digits_.begin() + static_cast<long>(digit_shift),
                          digits_.end());
  } else {
    const uint32_t mask = static_cast<uint32_t>((1ULL << bit_shift) - 1);
    uint32_t carry = 0;

    for (auto it = digits_.crbegin();
         it != digits_.crend() - static_cast<long>(digit_shift); ++it) {
      uint32_t shifted_digit = (*it >> bit_shift) | carry;
      carry = (*it & mask) << (32 - bit_shift);
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
  const bool this_zero = isZero();
  const bool rhs_zero = rhs.isZero();

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
    if (mag_cmp > 0) return std::strong_ordering::greater;
    if (mag_cmp < 0) return std::strong_ordering::less;
  } else {
    if (mag_cmp > 0) return std::strong_ordering::less;
    if (mag_cmp < 0) return std::strong_ordering::greater;
  }

  return std::strong_ordering::equal;
}

bool BigInt::operator==(const BigInt& rhs) const noexcept {
  return (*this <=> rhs) == std::strong_ordering::equal;
}

int BigInt::compareMagnitude(const BigInt& other) const noexcept {
  if (digits_.size() > other.digits_.size()) {
    return 1;
  } else if (digits_.size() < other.digits_.size()) {
    return -1;
  }

  for (auto left_it = digits_.crbegin(), right_it = other.digits_.crbegin();
       left_it != digits_.crend(); ++left_it, ++right_it) {
    if (*left_it > *right_it) return 1;
    if (*left_it < *right_it) return -1;
  }

  return 0;
}

BigInt abs(const BigInt& value) {
  BigInt result(value);
  result.positive_ = true;
  return result;
}

BigInt sqrt(const BigInt& value) {
  if (value.isZero()) {
    return BigInt();
  }
  if (!value.positive_) {
    throw std::domain_error("Square root of negative number is undefined");
  }
  if (value == BigInt(1)) {
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
  if (!a) {
    return b;
  } else if (!b) {
    return a;
  }

  BigInt result(1);
  a.positive_ = true;
  b.positive_ = true;

  while (a.isEven() && b.isEven()) {
    a >>= 1;
    b >>= 1;
    result <<= 1;
  }

  while (a.isEven()) {
    a >>= 1;
  }

  while (b) {
    while (b.isEven()) {
      b >>= 1;
    }
    BigInt temp = a;
    a = min(a, b);
    b = abs(temp - b);
  }

  result *= a;
  return result;
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
bool millerRabinWitness(const BigInt& witness,
                        const BigInt& d,
                        size_t s,
                        const BigInt& n) {
  const BigInt n_minus_1 = n - BigInt(1);
  BigInt x = powmod(witness, d, n);

  if (x == BigInt(1) || x == n_minus_1) {
    return true;
  }

  for (size_t r = 1; r < s; ++r) {
    x = powmod(x, BigInt(2), n);
    if (x == n_minus_1) {
      return true;
    }
    if (x == BigInt(1)) {
      return false;
    }
  }

  return false;
}

}  // anonymous namespace

bool BigInt::isProbablePrime(size_t rounds) const {
  if (!positive_ || isZero()) {
    return false;
  }
  if (*this == BigInt(2) || *this == BigInt(3)) {
    return true;
  }
  if (isEven() || *this < BigInt(2)) {
    return false;
  }

  // Quick divisibility check using small primes from anonymous namespace
  for (uint32_t p : kSmallPrimes) {
    if (*this == BigInt(p)) {
      return true;
    }
    if ((*this % BigInt(p)).isZero()) {
      return false;
    }
  }

  // Write n-1 as 2^s * d where d is odd
  BigInt d = *this - BigInt(1);
  size_t s = 0;
  while (d.isEven()) {
    d >>= 1;
    ++s;
  }

  // Deterministic witnesses for small numbers
  if (*this < BigInt(2047)) {
    return millerRabinWitness(BigInt(2), d, s, *this);
  }
  if (*this < BigInt(1373653)) {
    return millerRabinWitness(BigInt(2), d, s, *this) &&
           millerRabinWitness(BigInt(3), d, s, *this);
  }
  if (*this < BigInt(25326001)) {
    return millerRabinWitness(BigInt(2), d, s, *this) &&
           millerRabinWitness(BigInt(3), d, s, *this) &&
           millerRabinWitness(BigInt(5), d, s, *this);
  }
  if (*this < BigInt(3215031751ULL)) {
    return millerRabinWitness(BigInt(2), d, s, *this) &&
           millerRabinWitness(BigInt(3), d, s, *this) &&
           millerRabinWitness(BigInt(5), d, s, *this) &&
           millerRabinWitness(BigInt(7), d, s, *this);
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

  size_t num_words = (num_bits + 31) / 32;
  std::vector<uint32_t> digits(num_words);

  for (size_t i = 0; i < num_words; ++i) {
    digits[i] = dist32(gen);
  }

  size_t top_bits = num_bits % 32;
  if (top_bits == 0) {
    top_bits = 32;
  }

  uint32_t mask = (1U << top_bits) - 1;
  digits.back() &= mask;
  digits.back() |= (1U << (top_bits - 1));

  BigInt result;
  result.digits_ = std::move(digits);
  result.positive_ = true;
  return result;
}

BigInt BigInt::randomBelow(const BigInt& upper_bound) {
  if (upper_bound <= BigInt(1)) {
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

    if (candidate.isEven()) {
      candidate.digits_[0] |= 1;
    }

    bool divisible = false;
    for (uint32_t p : kSmallPrimes) {
      if (candidate == BigInt(p)) {
        return candidate;
      }
      if ((candidate % BigInt(p)).isZero()) {
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
  if (*this <= BigInt(2)) {
    return BigInt(2);
  }
  if (*this == BigInt(3)) {
    return BigInt(3);
  }

  BigInt candidate = *this;
  if (candidate.isEven()) {
    candidate += BigInt(1);
  }

  const size_t bits = candidate.bitLength();
  const size_t max_iterations =
      std::max(static_cast<size_t>(1000000), bits * bits * 100);
  size_t iterations = 0;

  while (iterations < max_iterations) {
    ++iterations;

    bool divisible = false;
    for (uint32_t p : kSmallPrimes) {
      if (candidate == BigInt(p)) {
        return candidate;
      }
      if ((candidate % BigInt(p)).isZero()) {
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

std::ostream& operator<<(std::ostream& out, const BigInt& value) {
  std::string str = value.toStdString(kDefaultOutputBase);
  out << str;
  return out;
}

std::istream& operator>>(std::istream& in, BigInt& value) {
  std::string str;
  in >> str;
  value = BigInt(str);
  return in;
}

BigInt BarrettReduction(const BigInt& dividend,
                        const BigInt& divisor,
                        const BigInt& mu) {
  BigInt remainder =
      dividend -
      ((dividend.shiftDigitsToLow(divisor.digits_.size() - 1) * mu)
           .shiftDigitsToLow(divisor.digits_.size() + 1) *
       divisor);

  while (remainder >= divisor) {
    remainder -= divisor;
  }

  return remainder;
}

std::string BigInt::toStdString(uint8_t base) const {
  std::stringstream ss;

  if (digits_.empty() || !(*this)) {
    return "0";
  }

  if (!positive_) {
    ss << '-';
  }

  if (base == kBaseBinary) {
    for (auto it = digits_.crbegin(); it != digits_.crend(); ++it) {
      ss << std::bitset<sizeof(uint32_t) * 8>(*it);
    }
  } else if (base == kBaseHexadecimal) {
    for (auto it = digits_.crbegin(); it != digits_.crend(); ++it) {
      ss << std::hex << std::setw(8) << std::setfill('0') << *it;
    }
  } else {  // base == kBaseDecimal
    const BigInt decimal_repr = toBigIntDec();
    for (auto it = decimal_repr.digits_.crbegin();
         it != decimal_repr.digits_.crend(); ++it) {
      ss << std::dec << std::setw(9) << std::setfill('0') << *it;
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

std::vector<uint32_t> BigInt::toStdVectorUint32_t() const {
  std::vector<uint32_t> result = digits_;
  std::reverse(result.begin(), result.end());
  return result;
}

std::vector<uint8_t> BigInt::toStdVectorUint8_t() const {
  std::vector<uint8_t> result;
  size_t num_bytes = byteLength();
  result.reserve(num_bytes);

  auto it = digits_.crbegin();
  size_t partial_bytes = num_bytes % sizeof(uint32_t);

  if (partial_bytes) {
    for (size_t i = 0; i < partial_bytes; ++i) {
      result.emplace_back(
          static_cast<uint8_t>(*it >> ((partial_bytes - i - 1) * 8)));
    }
    ++it;
  }

  while (it != digits_.crend()) {
    for (uint8_t i = 0; i < sizeof(uint32_t); ++i) {
      result.emplace_back(
          static_cast<uint8_t>(*it >> ((sizeof(uint32_t) - i - 1) * 8)));
    }
    ++it;
  }

  return result;
}

BigInt::operator uint64_t() const noexcept {
  if (digits_.empty()) {
    return 0;
  }
  if (digits_.size() >= 2) {
    return (static_cast<uint64_t>(*std::next(digits_.cbegin())) << 32) |
           static_cast<uint64_t>(*digits_.cbegin());
  }
  return digits_.front();
}

BigInt::operator uint32_t() const noexcept {
  return digits_.empty() ? 0 : digits_.front();
}

BigInt::operator uint16_t() const noexcept {
  return digits_.empty() ? static_cast<uint16_t>(0)
                         : static_cast<uint16_t>(digits_.front());
}

BigInt::operator uint8_t() const noexcept {
  return digits_.empty() ? static_cast<uint8_t>(0)
                         : static_cast<uint8_t>(digits_.front());
}

BigInt::operator int64_t() const noexcept {
  const auto magnitude = static_cast<uint64_t>(*this);
  return positive_ ? static_cast<int64_t>(magnitude)
                   : -static_cast<int64_t>(magnitude);
}

BigInt::operator int32_t() const noexcept {
  const auto magnitude = static_cast<uint32_t>(*this);
  return positive_ ? static_cast<int32_t>(magnitude)
                   : -static_cast<int32_t>(magnitude);
}

BigInt::operator int16_t() const noexcept {
  const auto magnitude = static_cast<uint16_t>(*this);
  return positive_ ? static_cast<int16_t>(magnitude)
                   : -static_cast<int16_t>(magnitude);
}

BigInt::operator int8_t() const noexcept {
  const auto magnitude = static_cast<uint8_t>(*this);
  return positive_ ? static_cast<int8_t>(magnitude)
                   : -static_cast<int8_t>(magnitude);
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

size_t BigInt::digitCount() const noexcept {
  return digits_.size();
}

bool BigInt::isZero() const noexcept {
  return digits_.empty() || (digits_.size() == 1 && digits_[0] == 0);
}

bool BigInt::isEven() const noexcept {
  return digits_.empty() || !(digits_.front() & 1);
}

bool BigInt::isOdd() const noexcept {
  return !digits_.empty() && (digits_.front() & 1);
}

bool BigInt::isPositive() const noexcept {
  return positive_ && !isZero();
}

bool BigInt::isNegative() const noexcept {
  return !positive_ && !isZero();
}

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
    result.digits_.erase(result.digits_.begin(),
                         result.digits_.begin() + static_cast<long>(shift));
  } else {
    result.digits_.shrink_to_fit();
    result.digits_.reserve(1);
    result.digits_.emplace_back(0);
    result.positive_ = true;
  }

  return result;
}

BigInt BigInt::multiplySchoolbook(const BigInt& other) const {
  if (digits_.empty() || other.digits_.empty() || isZero() || other.isZero()) {
    return BigInt();
  }

  const size_t m = digits_.size();
  const size_t n = other.digits_.size();

  BigInt result;
  result.digits_.resize(m + n, 0);

  for (size_t i = 0; i < m; ++i) {
    uint64_t carry = 0;
    for (size_t j = 0; j < n; ++j) {
      uint64_t product = static_cast<uint64_t>(digits_[i]) *
                             static_cast<uint64_t>(other.digits_[j]) +
                             static_cast<uint64_t>(result.digits_[i + j]) +
                         carry;
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
  if (isZero() || other.isZero()) {
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
    low1.digits_.assign(digits_.begin(),
                        digits_.begin() + static_cast<long>(half));
    low1.positive_ = true;
    low1.deleteZeroHighOrderDigit();
    high1.digits_.assign(digits_.begin() + static_cast<long>(half),
                         digits_.end());
    high1.positive_ = true;
    high1.deleteZeroHighOrderDigit();
  }

  // Split other = high2 * B^half + low2
  BigInt low2, high2;
  if (n <= half) {
    low2 = other;
    low2.positive_ = true;
  } else {
    low2.digits_.assign(other.digits_.begin(),
                        other.digits_.begin() + static_cast<long>(half));
    low2.positive_ = true;
    low2.deleteZeroHighOrderDigit();
    high2.digits_.assign(other.digits_.begin() + static_cast<long>(half),
                         other.digits_.end());
    high2.positive_ = true;
    high2.deleteZeroHighOrderDigit();
  }

  // Karatsuba: 3 multiplications instead of 4
  BigInt z0 = low1.multiplyKaratsuba(low2);
  BigInt z2 = high1.multiplyKaratsuba(high2);

  BigInt sum1 = abs(low1) + abs(high1);
  BigInt sum2 = abs(low2) + abs(high2);
  BigInt z1 = sum1.multiplyKaratsuba(sum2) - z0 - z2;

  BigInt result =
      z2.shiftDigitsToHigh(2 * half) + z1.shiftDigitsToHigh(half) + z0;
  result.positive_ = true;
  return result;
}

BigInt BigInt::toBigIntDec() const {
  const BigInt kDecimalBase(1000000000);
  BigInt value = abs(*this);
  BigInt result;
  result.positive_ = positive_;
  result.digits_.reserve(digits_.size() + 1);

  while (value) {
    auto [quotient, remainder] = value.DivMod(kDecimalBase);
    result.digits_.emplace_back(remainder.digits_.front());
    value = quotient;
  }

  return result;
}

}  // namespace bigint



