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

BigInt BigInt::operator +() const
{
    return *this;
}

BigInt BigInt::operator + (const BigInt& addend) const
{
    if(positive_ && addend.positive_)
    {
        BigInt sum;
        uint32_t carry = 0;
        uint64_t sum_temp;
        bool augendGreater = (digits_.size() >= addend.digits_.size());
        sum.digits_.reserve(augendGreater ? digits_.size() + 1 : addend.digits_.size() + 1);
        std::vector<uint32_t>::const_iterator iteratorAugend = augendGreater ? digits_.cbegin() : addend.digits_.cbegin();
        std::vector<uint32_t>::const_iterator iteratorAddend = augendGreater ? addend.digits_.cbegin() : digits_.cbegin();
        std::vector<uint32_t>::const_iterator iteratorAugendEnd = augendGreater ? digits_.cend() : addend.digits_.cend();
        std::vector<uint32_t>::const_iterator iteratorAddendEnd = augendGreater ? addend.digits_.cend() : digits_.cend();
        while(iteratorAddend != iteratorAddendEnd)
        {
            sum_temp = static_cast<uint64_t>(*iteratorAugend) + static_cast<uint64_t>(*iteratorAddend) + static_cast<uint64_t>(carry);
            sum.digits_.emplace_back(static_cast<uint32_t>(sum_temp & UINT32_MAX));
            carry = static_cast<uint32_t>(sum_temp >> 32);
            ++iteratorAugend;
            ++iteratorAddend;
        }
        while(iteratorAugend != iteratorAugendEnd)
        {
            sum_temp = static_cast<uint64_t>(*iteratorAugend) + static_cast<uint64_t>(carry);
            sum.digits_.emplace_back(static_cast<uint32_t>(sum_temp & UINT32_MAX));
            carry = static_cast<uint32_t>(sum_temp >> 32);
            ++iteratorAugend;
        }
        if(carry)
        {
            sum.digits_.emplace_back(carry);
        }
        sum.positive_ = true;
        return sum;
    }
    else if(positive_ && !addend.positive_)
    {
        return *this - abs(addend);
    }
    else if(!positive_ && addend.positive_)
    {
        return addend - abs(*this);
    }
    else // !positive_ && !addend.positive_
    {
        return -(abs(*this) + abs(addend));
    }
}

BigInt& BigInt::operator += (const BigInt& augend)
{
    // Handle zero cases
    if (augend.isZero()) {
        return *this;
    }
    if (isZero()) {
        *this = augend;
        return *this;
    }

    // Same sign: add magnitudes
    if (positive_ == augend.positive_) {
        // Ensure we have enough space
        if (digits_.size() < augend.digits_.size()) {
            digits_.resize(augend.digits_.size(), 0);
        }

        uint32_t carry = 0;
        size_t i = 0;
        for (; i < augend.digits_.size(); ++i) {
            uint64_t sum = static_cast<uint64_t>(digits_[i]) +
                           static_cast<uint64_t>(augend.digits_[i]) +
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
        // Compare magnitudes
        int cmp = compareMagnitude(augend);
        if (cmp == 0) {
            // Equal magnitudes, result is zero
            digits_.clear();
            digits_.push_back(0);
            positive_ = true;
            return *this;
        }

        const BigInt* larger;
        const BigInt* smaller;
        bool resultPositive;

        if (cmp > 0) {
            // |*this| > |augend|
            larger = this;
            smaller = &augend;
            resultPositive = positive_;
        } else {
            // |*this| < |augend|
            larger = &augend;
            smaller = this;
            resultPositive = augend.positive_;
        }

        // Subtract smaller from larger
        std::vector<uint32_t> result;
        result.resize(larger->digits_.size(), 0);

        uint8_t borrow = 0;
        for (size_t i = 0; i < larger->digits_.size(); ++i) {
            int64_t diff = static_cast<int64_t>(larger->digits_[i]) -
                           (i < smaller->digits_.size() ? static_cast<int64_t>(smaller->digits_[i]) : 0) -
                           borrow;
            if (diff >= 0) {
                result[i] = static_cast<uint32_t>(diff);
                borrow = 0;
            } else {
                result[i] = static_cast<uint32_t>(diff + static_cast<int64_t>(kBasisCalcSys));
                borrow = 1;
            }
        }

        digits_ = std::move(result);
        positive_ = resultPositive;
        deleteZeroHighOrderDigit();
    }

    return *this;
}

BigInt& BigInt::operator ++()
{
    if (positive_) {
        // Positive number: add 1
        if (digits_.empty()) {
            digits_.push_back(1);
            return *this;
        }
        for (size_t i = 0; i < digits_.size(); ++i) {
            if (digits_[i] < UINT32_MAX) {
                ++digits_[i];
                return *this;  // No carry, done
            }
            digits_[i] = 0;  // Carry to next digit
        }
        digits_.push_back(1);  // All digits overflowed, add new digit
    } else {
        // Negative number: subtract 1 from magnitude
        for (size_t i = 0; i < digits_.size(); ++i) {
            if (digits_[i] > 0) {
                --digits_[i];
                deleteZeroHighOrderDigit();
                if (isZero()) positive_ = true;  // Normalize -0 to +0
                return *this;
            }
            digits_[i] = UINT32_MAX;  // Borrow from next digit
        }
    }
    return *this;
}

BigInt BigInt::operator ++(int)
{
    const BigInt bigNum = *this;
    ++(*this);
    return bigNum;
}

BigInt BigInt::operator -() const
{
    BigInt negative = *this;
    negative.positive_ = !positive_;
    return negative;
}

BigInt BigInt::operator - (const BigInt& subtrahend) const
{
    if(positive_ && subtrahend.positive_)
    {
        if(*this >= subtrahend)
        {
            BigInt difference;
            uint8_t borrow = 0;
            int64_t difference_temp;
            difference.digits_.reserve(digits_.size());
            std::vector<uint32_t>::const_iterator iteratorMinuend = digits_.cbegin();
            std::vector<uint32_t>::const_iterator iteratorSubtrahend = subtrahend.digits_.cbegin();
            while(iteratorSubtrahend != subtrahend.digits_.cend())
            {
                difference_temp = static_cast<int64_t>(*iteratorMinuend) - static_cast<int64_t>(*iteratorSubtrahend) - static_cast<int64_t>(borrow);
                if(difference_temp >= 0)
                {
                    difference.digits_.emplace_back(static_cast<uint32_t>(difference_temp));
                    borrow = 0;
                }
                else // difference_temp < 0
                {
                    difference.digits_.emplace_back(static_cast<uint32_t>(difference_temp + static_cast<int64_t>(kBasisCalcSys)));
                    borrow = 1;
                }
                ++iteratorMinuend;
                ++iteratorSubtrahend;
            }
            while(iteratorMinuend != digits_.cend())
            {
                difference_temp = static_cast<int64_t>(*iteratorMinuend) - static_cast<int64_t>(borrow);
                if(difference_temp >= 0)
                {
                    difference.digits_.emplace_back(static_cast<uint32_t>(difference_temp));
                    borrow = 0;
                }
                else // difference_temp < 0
                {
                    difference.digits_.emplace_back(static_cast<uint32_t>(difference_temp + static_cast<int64_t>(kBasisCalcSys)));
                    borrow = 1;
                }
                ++iteratorMinuend;
            }
            difference.deleteZeroHighOrderDigit();
            difference.positive_ = true;
            return difference;
        }
        else // minuend < subtrahend
        {
            return -(subtrahend - *this);
        }
    }
    else if(!positive_ && subtrahend.positive_)
    {
        return -(abs(*this) + subtrahend);
    }
    else if(positive_ && !subtrahend.positive_)
    {
        return *this + abs(subtrahend);
    }
    else // !positive_ && !subtrahend.positive_
    {
        return abs(subtrahend) - abs(*this);
    }
}

BigInt& BigInt::operator -= (const BigInt& subtrahend)
{
    // a -= b is equivalent to a += (-b)
    // Create a temporary with flipped sign and use +=
    if (subtrahend.isZero()) {
        return *this;
    }

    BigInt negated = subtrahend;
    negated.positive_ = !subtrahend.positive_;
    return *this += negated;
}

BigInt& BigInt::operator -- ()
{
    if (positive_) {
        // Positive number: subtract 1 from magnitude
        if (digits_.empty() || isZero()) {
            // 0 - 1 = -1
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
            digits_[i] = UINT32_MAX;  // Borrow from next digit
        }
    } else {
        // Negative number: add 1 to magnitude
        if (digits_.empty()) {
            digits_.push_back(1);
            return *this;
        }
        for (size_t i = 0; i < digits_.size(); ++i) {
            if (digits_[i] < UINT32_MAX) {
                ++digits_[i];
                return *this;  // No carry, done
            }
            digits_[i] = 0;  // Carry to next digit
        }
        digits_.push_back(1);  // All digits overflowed, add new digit
    }
    return *this;
}

BigInt BigInt::operator -- (int)
{
    const BigInt bigNum = *this;
    --(*this);
    return bigNum;
}

BigInt BigInt::operator * (const uint32_t multiplier) const
{
    BigInt product;
    uint32_t carry = 0;
    uint64_t product_temp;
    product.digits_.reserve(digits_.size() + 1);
    for(std::vector<uint32_t>::const_iterator iteratorMultiplicand = digits_.cbegin(); iteratorMultiplicand != digits_.cend(); ++iteratorMultiplicand)
    {
        product_temp = static_cast<uint64_t>(*iteratorMultiplicand) * static_cast<uint64_t>(multiplier) + static_cast<uint64_t>(carry);
        product.digits_.emplace_back(static_cast<uint32_t>(product_temp & UINT32_MAX));
        carry = static_cast<uint32_t>(product_temp >> 32);
    }
    if(carry)
    {
        product.digits_.emplace_back(carry);
    }
    product.positive_ = positive_;
    return product;
}

BigInt& BigInt::operator *= (const uint32_t multiplier)
{
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

BigInt BigInt::operator * (const BigInt& multiplier) const
{
    // Handle zero cases early
    if (isZero() || multiplier.isZero()) {
        return BigInt();
    }

    // Use Karatsuba for large numbers, schoolbook for small
    BigInt product;
    const size_t maxSize = std::max(digits_.size(), multiplier.digits_.size());
    if (maxSize >= kKaratsubaThreshold) {
        product = multiplyKaratsuba(multiplier);
    } else {
        product = multiplySchoolbook(multiplier);
    }

    // Set the sign: positive if both have same sign
    product.positive_ = (positive_ == multiplier.positive_);
    return product;
}

BigInt& BigInt::operator *= (const BigInt& multiplier)
{
    *this = *this * multiplier;
    return *this;
}

std::pair<BigInt, BigInt> BigInt::DivMod(const BigInt& divisor) const
{
    if (divisor.isZero()) {
        throw std::domain_error("Division by zero");
    }

    const size_t bitLengthDivisor = divisor.bitLength();
    BigInt fraction(0);
    BigInt remainder(abs(*this));
    BigInt borrow;
    size_t differenceRemainderNDivisorbitLength;
    fraction.digits_.reserve(digits_.size());
    while(remainder >= abs(divisor))
    {
        differenceRemainderNDivisorbitLength = remainder.bitLength() - bitLengthDivisor;
        borrow = abs(divisor) << differenceRemainderNDivisorbitLength;
        if(remainder < borrow)
        {
            borrow >>= 1;
            --differenceRemainderNDivisorbitLength;
        }
        remainder -= borrow;
        fraction += BigInt(1) << differenceRemainderNDivisorbitLength;  // 1 << n = 2^n
    }
    fraction.positive_ = positive_ == divisor.positive_;
    remainder.positive_ = positive_;
    if(divisor.positive_)
    {
        while(!remainder.positive_)
        {
            remainder += divisor;
        }
    }
    return std::make_pair(fraction, remainder);
}

BigInt BigInt::operator / (const BigInt& divisor) const
{
    return DivMod(divisor).first;
}

BigInt& BigInt::operator /= (const BigInt& divisor)
{
    *this = DivMod(divisor).first;
    return *this;
}

BigInt BigInt::operator % (const BigInt& divisor) const
{
    return DivMod(divisor).second;
}

BigInt& BigInt::operator %= (const BigInt& divisor)
{
    *this = DivMod(divisor).second;
    return *this;
}

BigInt pow(const BigInt& base, const BigInt& exponent)
{
    if(!exponent.positive_)
    {
        return BigInt(0);
    }
    BigInt power(1);
    power.digits_.reserve(base.digits_.size() * static_cast<size_t>(exponent));
    for(size_t indexBitExponent = exponent.bitLength() - 1; indexBitExponent > 0; --indexBitExponent)
    {
        if(exponent.digits_[indexBitExponent >> 5] & (1 << (indexBitExponent & 31)))
        {
            power *= base;
        }
        power *= power;
    }
    if(exponent.isOdd())
    {
        power *= base;
    }
    if(!base.positive_)
    {
        power.positive_ = exponent.isEven();
    }
    return power;
}

size_t log2(const BigInt& antilogarithm) noexcept
{
    return antilogarithm.bitLength() - 1;
}

/*BigInt powmod(BigInt base, const BigInt& exponent, const BigInt& divisor)
{
    BigInt power(1);
    power.digits_.reserve(divisor.digits_.size());
    const BigInt mu = power.shiftDigitsToHigh(divisor.digits_.size() * 2) / divisor;
    const uint32_t bitLengthExponent = exponent.bitLength();
    for(size_t indexBitExponent = 0; indexBitExponent < bitLengthExponent; ++indexBitExponent)
    {
        if(exponent.digits_[indexBitExponent >> 5] & (1 << (indexBitExponent & 31)))
        {
            power = BarrettReduction(power * base, divisor, mu);
        }
        base = BarrettReduction(base * base, divisor, mu);
    }
    return power;
}*/

BigInt powmod(const BigInt& base, const BigInt& exponent, const BigInt& divisor)
{
    if (divisor.isZero()) {
        throw std::domain_error("Modular exponentiation with zero divisor");
    }

    BigInt power(1);
    power.digits_.reserve(divisor.digits_.size());
    const size_t bitLen = exponent.bitLength();
    for(size_t i = 0; i < bitLen; ++i)
    {
        size_t indexBitExponent = bitLen - 1 - i;
        power = (power * power) % divisor;
        if(exponent.digits_[indexBitExponent >> 5] & (1 << (indexBitExponent & 31)))
        {
            power = (power * base) % divisor;
        }
    }
    return power;
}

BigInt inversemod(BigInt dividend, const BigInt& divisor)
{
    if (divisor.isZero()) {
        throw std::domain_error("Modular inverse with zero divisor");
    }
    if (!isCoprime(dividend, divisor)) {
        throw std::domain_error("Modular inverse does not exist (numbers not coprime)");
    }

    BigInt divisor_copy(divisor);
    BigInt fraction;
    BigInt x0(0);
    BigInt x1(1);
    BigInt x_temp;
    while(dividend > BigInt(1))
    {
        fraction = dividend / divisor_copy;
        x_temp = divisor_copy;
        divisor_copy = dividend % divisor_copy;
        dividend = x_temp;
        x_temp = x0;
        x0 = x1 - (fraction * x0);
        x1 = x_temp;
    }
    if(!x1.positive_)
    {
        x1 += divisor;
    }
    return x1;
}

bool congruencemod(const BigInt& dividend1, const BigInt& dividend2, const BigInt& divisor)
{
    if (divisor.isZero()) {
        throw std::domain_error("Congruence check with zero divisor");
    }

    BigInt remainder1(dividend1 % divisor);
    BigInt remainder2(dividend2 % divisor);
    while(!remainder1.positive_)
    {
        remainder1 += divisor;
    }
    while(remainder1 > divisor)
    {
        remainder1 -= divisor;
    }
    while(!remainder2.positive_)
    {
        remainder2 += divisor;
    }
    while(remainder2 > divisor)
    {
        remainder2 -= divisor;
    }
    return remainder1 == remainder2;
}

bool isCoprime(const BigInt& bigInt1, const BigInt& bigInt2)
{
    return gcd(bigInt1, bigInt2) == BigInt(1);
}

int8_t symbolJacobi(BigInt bigInt1, BigInt bigInt2)
{
    if(!isCoprime(bigInt1, bigInt2))
    {
        return 0;
    }
    int8_t result = 1;
    size_t iterator = 0;
    BigInt bigInt3;
    if(!bigInt1.positive_)
    {
        bigInt1.positive_ = true;
        if(bigInt2 % BigInt(4) == BigInt(3))
        {
            result = -result;
        }
    }
    while(bigInt1)
    {
        iterator = 0;
        while(bigInt1.isEven())
        {
            bigInt1 >>= 1;
            ++iterator;
        }
        if(iterator % 2)
        {
            if(bigInt2 % BigInt(8) == BigInt(3) || bigInt2 % BigInt(8) == BigInt(5))
            {
                result = -result;
            }
        }
        if(bigInt1 % BigInt(4) == BigInt(3) && bigInt2 % BigInt(4) == BigInt(3))
        {
            result = -result;
        }
        bigInt3 = bigInt1;
        bigInt1 = bigInt2 % bigInt3;
        bigInt2 = bigInt3;
    }
    return result;
}

BigInt BigInt::operator ~() const
{
    return -*this - BigInt(1);
}

BigInt BigInt::operator & (const BigInt& rightBitwiseAND) const
{
    if(positive_ && rightBitwiseAND.positive_)
    {
        BigInt bitwiseAND;
        bitwiseAND.digits_.reserve(std::min(digits_.size(), rightBitwiseAND.digits_.size()));
        std::vector<uint32_t>::const_iterator iteratorLeftBitwiseAND = digits_.cbegin();
        std::vector<uint32_t>::const_iterator iteratorRightBitwiseAND = rightBitwiseAND.digits_.cbegin();
        while(iteratorLeftBitwiseAND != digits_.cend() && iteratorRightBitwiseAND != rightBitwiseAND.digits_.cend())
        {
            bitwiseAND.digits_.emplace_back(*iteratorLeftBitwiseAND & *iteratorRightBitwiseAND);
            ++iteratorLeftBitwiseAND;
            ++iteratorRightBitwiseAND;
        }
        bitwiseAND.positive_ = true;
        return bitwiseAND;
    }
    else if(!positive_ && !rightBitwiseAND.positive_)
    {
        return -(~*this | ~rightBitwiseAND) - BigInt(1);
    }
    else if(positive_ && !rightBitwiseAND.positive_)
    {
        return (*this | ~rightBitwiseAND) + rightBitwiseAND + BigInt(1);
    }
    else // !positive_ && rightBitwiseAnd.positive_
    {
        return (~*this | rightBitwiseAND) + *this + BigInt(1);
    }
}

BigInt& BigInt::operator &= (const BigInt& rightBitwiseAND)
{
    *this = *this & rightBitwiseAND;
    return *this;
}

BigInt BigInt::operator | (const BigInt& rightBitwiseOR) const
{
    if(positive_ && rightBitwiseOR.positive_)
    {
        BigInt bitwiseOR;
        const bool leftBitwiseORGreater = (digits_.size() >= rightBitwiseOR.digits_.size());
        bitwiseOR.digits_.reserve(leftBitwiseORGreater ? digits_.size() + 1 : rightBitwiseOR.digits_.size() + 1);
        std::vector<uint32_t>::const_iterator iteratorLeftBitwiseOR = leftBitwiseORGreater ? digits_.cbegin() : rightBitwiseOR.digits_.cbegin();
        std::vector<uint32_t>::const_iterator iteratorRightBitwiseOR = leftBitwiseORGreater ? rightBitwiseOR.digits_.cbegin() : digits_.cbegin();
        std::vector<uint32_t>::const_iterator iteratorLeftBitwiseOREnd = leftBitwiseORGreater ? digits_.cend() : rightBitwiseOR.digits_.cend();
        std::vector<uint32_t>::const_iterator iteratorRightBitwiseOREnd = leftBitwiseORGreater ? rightBitwiseOR.digits_.cend() : digits_.cend();
        while(iteratorRightBitwiseOR != iteratorRightBitwiseOREnd)
        {
            bitwiseOR.digits_.emplace_back(*iteratorLeftBitwiseOR | *iteratorRightBitwiseOR);
            ++iteratorLeftBitwiseOR;
            ++iteratorRightBitwiseOR;
        }
        while(iteratorLeftBitwiseOR != iteratorLeftBitwiseOREnd)
        {
            bitwiseOR.digits_.emplace_back(*iteratorLeftBitwiseOR);
            ++iteratorLeftBitwiseOR;
        }
        bitwiseOR.positive_ = true;
        return bitwiseOR;
    }
    else if(!positive_ && !rightBitwiseOR.positive_)
    {
        return -(~*this & ~rightBitwiseOR) - BigInt(1);
    }
    else if(positive_ && !rightBitwiseOR.positive_)
    {
        return (*this & ~rightBitwiseOR) + rightBitwiseOR;
    }
    else // !positive_ && rightBitwiseOR.positive_
    {
        return (~*this & rightBitwiseOR) + *this;
    }
}

BigInt& BigInt::operator |= (const BigInt& rightBitwiseOR)
{
    *this = *this | rightBitwiseOR;
    return *this;
}

BigInt BigInt::operator ^ (const BigInt& rightBitwiseXOR) const
{
    if(positive_ && rightBitwiseXOR.positive_)
    {
        BigInt bitwiseXOR;
        const bool leftBitwiseXORGreater = (digits_.size() >= rightBitwiseXOR.digits_.size());
        bitwiseXOR.digits_.reserve(leftBitwiseXORGreater ? digits_.size() + 1 : rightBitwiseXOR.digits_.size() + 1);
        std::vector<uint32_t>::const_iterator iteratorLeftBitwiseXOR = leftBitwiseXORGreater ? digits_.cbegin() : rightBitwiseXOR.digits_.cbegin();
        std::vector<uint32_t>::const_iterator iteratorRightBitwiseXOR = leftBitwiseXORGreater ? rightBitwiseXOR.digits_.cbegin() : digits_.cbegin();
        std::vector<uint32_t>::const_iterator iteratorLeftBitwiseXOREnd = leftBitwiseXORGreater ? digits_.cend() : rightBitwiseXOR.digits_.cend();
        std::vector<uint32_t>::const_iterator iteratorRightBitwiseXOREnd = leftBitwiseXORGreater ? rightBitwiseXOR.digits_.cend() : digits_.cend();
        while(iteratorRightBitwiseXOR != iteratorRightBitwiseXOREnd)
        {
            bitwiseXOR.digits_.emplace_back(*iteratorLeftBitwiseXOR ^ *iteratorRightBitwiseXOR);
            ++iteratorLeftBitwiseXOR;
            ++iteratorRightBitwiseXOR;
        }
        while(iteratorLeftBitwiseXOR != iteratorLeftBitwiseXOREnd)
        {
            bitwiseXOR.digits_.emplace_back(*iteratorLeftBitwiseXOR);
            ++iteratorLeftBitwiseXOR;
        }
        bitwiseXOR.positive_ = true;
        return bitwiseXOR;
    }
    else // !positive_ || !rightBitwiseXOR.positive_
    {
        return (*this | rightBitwiseXOR) & (~*this | ~rightBitwiseXOR);
    }
}

BigInt& BigInt::operator ^= (const BigInt& rightBitwiseXOR)
{
    *this = *this ^ rightBitwiseXOR;
    return *this;
}

BigInt BigInt::operator << (const size_t shift) const
{
    // Handle zero shift or zero value early
    if(!shift || digits_.empty() || !(*this))
    {
        return *this;
    }

    // Split shift into digit shift (multiples of 32) and bit shift (remainder)
    const size_t digitShift = shift / 32;
    const size_t bitShift = shift % 32;

    BigInt shifted;
    shifted.positive_ = positive_;
    shifted.digits_.reserve(digits_.size() + digitShift + 1);

    // Insert zeros for digit shift (these become the low-order digits)
    shifted.digits_.resize(digitShift, 0);

    // Perform bit shift on the remaining bits
    if(bitShift == 0)
    {
        // No bit shift needed, just copy digits
        shifted.digits_.insert(shifted.digits_.end(), digits_.begin(), digits_.end());
    }
    else
    {
        uint32_t carry = 0;
        for(auto it = digits_.cbegin(); it != digits_.cend(); ++it)
        {
            uint32_t shifted_temp = (*it << bitShift) | carry;
            carry = *it >> (32 - bitShift);
            shifted.digits_.emplace_back(shifted_temp);
        }
        if(carry)
        {
            shifted.digits_.emplace_back(carry);
        }
    }
    return shifted;
}

BigInt& BigInt::operator <<= (const size_t shift)
{
    *this = *this << shift;
    return *this;
}

BigInt BigInt::operator >> (const size_t shift) const
{
    // Handle zero shift or zero value early
    if(!shift || digits_.empty() || !(*this))
    {
        return *this;
    }

    // Split shift into digit shift (multiples of 32) and bit shift (remainder)
    const size_t digitShift = shift / 32;
    const size_t bitShift = shift % 32;

    // If we're shifting away all digits, return zero
    if(digitShift >= digits_.size())
    {
        return BigInt(0);
    }

    BigInt shifted;
    shifted.positive_ = positive_;

    // Skip the low-order digits that are shifted away
    const size_t remainingDigits = digits_.size() - digitShift;
    shifted.digits_.reserve(remainingDigits);

    if(bitShift == 0)
    {
        // No bit shift needed, just copy remaining digits
        shifted.digits_.insert(shifted.digits_.end(),
                               digits_.begin() + static_cast<long>(digitShift),
                               digits_.end());
    }
    else
    {
        const uint32_t mask = static_cast<uint32_t>((1ULL << bitShift) - 1);
        // Process from high to low, building result in reverse order
        uint32_t carry = 0;
        for(auto it = digits_.crbegin(); it != digits_.crend() - static_cast<long>(digitShift); ++it)
        {
            uint32_t shifted_temp = (*it >> bitShift) | carry;
            carry = (*it & mask) << (32 - bitShift);
            shifted.digits_.emplace_back(shifted_temp);
        }
        // Reverse to get correct order (low to high)
        std::reverse(shifted.digits_.begin(), shifted.digits_.end());
    }

    shifted.deleteZeroHighOrderDigit();
    return shifted;
}

BigInt& BigInt::operator >>= (const size_t shift)
{
    *this = *this >> shift;
    return *this;
}

BigInt BigInt::leftCircularShift(const size_t shift) const
{
    const BigInt mask(--(BigInt(1) << bitLength()));
    return (((*this << shift) | (*this >> (bitLength() - shift))) & mask);
}

BigInt BigInt::rightCircularShift(const size_t shift) const
{
    const BigInt mask(--(BigInt(1) << bitLength()));
    return (((*this >> shift) | (*this << (bitLength() - shift))) & mask);
}

bool BigInt::operator !() const noexcept
{
    for(const auto& digit : digits_)
    {
        if(digit != 0)
        {
            return false;
        }
    }
    return true;
}

bool BigInt::operator && (const BigInt& rightAND) const noexcept
{
    return static_cast<bool>(*this) && static_cast<bool>(rightAND);
}

bool BigInt::operator || (const BigInt& rightOR) const noexcept
{
    return static_cast<bool>(*this) || static_cast<bool>(rightOR);
}

std::strong_ordering BigInt::operator<=>(const BigInt& rhs) const noexcept
{
    // Handle zero cases: -0 == +0
    const bool leftZero = isZero();
    const bool rightZero = rhs.isZero();
    if (leftZero && rightZero) return std::strong_ordering::equal;

    // Different signs
    if (positive_ && !rhs.positive_) return std::strong_ordering::greater;
    if (!positive_ && rhs.positive_) return std::strong_ordering::less;

    // Same sign - compare magnitudes
    const int cmp = compareMagnitude(rhs);
    if (positive_)
    {
        // Both positive: larger magnitude = greater
        if (cmp > 0) return std::strong_ordering::greater;
        if (cmp < 0) return std::strong_ordering::less;
    }
    else
    {
        // Both negative: larger magnitude = less (e.g., -5 < -3)
        if (cmp > 0) return std::strong_ordering::less;
        if (cmp < 0) return std::strong_ordering::greater;
    }
    return std::strong_ordering::equal;
}

bool BigInt::operator==(const BigInt& rhs) const noexcept
{
    return (*this <=> rhs) == std::strong_ordering::equal;
}

int BigInt::compareMagnitude(const BigInt& other) const noexcept
{
    if (digits_.size() > other.digits_.size())
    {
        return 1;
    }
    else if (digits_.size() < other.digits_.size())
    {
        return -1;
    }
    // Same size - compare digit by digit from high to low
    for (auto itLeft = digits_.crbegin(), itRight = other.digits_.crbegin();
         itLeft != digits_.crend(); ++itLeft, ++itRight)
    {
        if (*itLeft > *itRight) return 1;
        if (*itLeft < *itRight) return -1;
    }
    return 0;
}

BigInt abs(const BigInt& bigInt)
{
    BigInt absolute(bigInt);
    absolute.positive_ = true;
    return absolute;
}

BigInt sqrt(const BigInt& value)
{
    // Handle special cases
    if (value.isZero()) return BigInt();
    if (!value.positive_) {
        throw std::domain_error("Square root of negative number is undefined");
    }
    if (value == BigInt(1)) return BigInt(1);

    // Newton's method (Heron's method) for integer square root
    // Initial guess: 2^((bitLength+1)/2) which is close to sqrt(value)
    size_t bits = value.bitLength();
    BigInt x = BigInt(1) << ((bits + 1) / 2);

    // Iterate until convergence
    BigInt x_prev;
    while (true) {
        // x_new = (x + value/x) / 2
        BigInt x_new = (x + value / x) >> 1;

        // Newton's method for integer sqrt converges when x_new >= x
        if (x_new >= x) {
            break;
        }
        x = x_new;
    }

    return x;
}

BigInt gcd(BigInt bigInt1, BigInt bigInt2)
{
    if(!bigInt1)
    {
        return bigInt2;
    }
    else if(!bigInt2)
    {
        return bigInt1;
    }
    BigInt greatestCommonDivisor(1);
    BigInt bigInt1_temp;
    bigInt1.positive_ = true;
    bigInt2.positive_ = true;
    while(bigInt1.isEven() && bigInt2.isEven())
    {
        bigInt1 >>= 1;
        bigInt2 >>= 1;
        greatestCommonDivisor <<= 1;
    }
    while(bigInt1.isEven())
    {
        bigInt1 >>= 1;
    }
    while(bigInt2)
    {
        while(bigInt2.isEven())
        {
            bigInt2 >>= 1;
        }
        bigInt1_temp = bigInt1;
        bigInt1 = min(bigInt1, bigInt2);
        bigInt2 = abs(bigInt1_temp - bigInt2);
    }
    greatestCommonDivisor *= bigInt1;
    return greatestCommonDivisor;
}

BigInt lcm(const BigInt& bigInt1, const BigInt& bigInt2)
{
    return (bigInt1 * bigInt2) / gcd(bigInt1, bigInt2);
}

const BigInt& max(const BigInt& bigInt1, const BigInt& bigInt2) noexcept
{
    return bigInt1 > bigInt2 ? bigInt1 : bigInt2;
}

const BigInt& min(const BigInt& bigInt1, const BigInt& bigInt2) noexcept
{
    return bigInt1 < bigInt2 ? bigInt1 : bigInt2;
}

void swap(BigInt& lhs, BigInt& rhs) noexcept
{
    std::swap(lhs.positive_, rhs.positive_);
    std::swap(lhs.digits_, rhs.digits_);
}

// ============================================================================
// Number Theory
// ============================================================================

namespace {

// Miller-Rabin witness test helper
bool millerRabinWitness(const BigInt& a, const BigInt& d, size_t s, const BigInt& n) {
    const BigInt n_minus_1 = n - BigInt(1);

    // Compute x = a^d mod n
    BigInt x = powmod(a, d, n);

    if (x == BigInt(1) || x == n_minus_1) {
        return true;  // Passes for this witness
    }

    for (size_t r = 1; r < s; ++r) {
        x = powmod(x, BigInt(2), n);
        if (x == n_minus_1) {
            return true;  // Passes
        }
        if (x == BigInt(1)) {
            return false;  // Composite
        }
    }

    return false;  // Composite
}

}  // anonymous namespace

bool BigInt::isProbablePrime(size_t rounds) const
{
    // Handle small cases
    if (!positive_ || isZero()) return false;
    if (*this == BigInt(2)) return true;
    if (*this == BigInt(3)) return true;
    if (isEven()) return false;
    if (*this < BigInt(2)) return false;

    // Small primes for quick divisibility check
    static const uint32_t smallPrimes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    for (uint32_t p : smallPrimes) {
        if (*this == BigInt(p)) return true;
        BigInt remainder = *this % BigInt(p);
        if (remainder.isZero()) return false;
    }

    // Write n-1 as 2^s * d where d is odd
    BigInt d = *this - BigInt(1);
    size_t s = 0;
    while (d.isEven()) {
        d >>= 1;
        ++s;
    }

    // Use deterministic witnesses for small numbers
    // For n < 2047, witness {2} is sufficient
    if (*this < BigInt(2047)) {
        return millerRabinWitness(BigInt(2), d, s, *this);
    }
    // For n < 1373653, witnesses {2, 3} are sufficient
    if (*this < BigInt(1373653)) {
        return millerRabinWitness(BigInt(2), d, s, *this) &&
               millerRabinWitness(BigInt(3), d, s, *this);
    }
    // For n < 25326001, witnesses {2, 3, 5} are sufficient
    if (*this < BigInt(25326001)) {
        return millerRabinWitness(BigInt(2), d, s, *this) &&
               millerRabinWitness(BigInt(3), d, s, *this) &&
               millerRabinWitness(BigInt(5), d, s, *this);
    }
    // For n < 3215031751, witnesses {2, 3, 5, 7} are sufficient
    if (*this < BigInt(3215031751ULL)) {
        return millerRabinWitness(BigInt(2), d, s, *this) &&
               millerRabinWitness(BigInt(3), d, s, *this) &&
               millerRabinWitness(BigInt(5), d, s, *this) &&
               millerRabinWitness(BigInt(7), d, s, *this);
    }

    // For larger numbers, use random witnesses (probabilistic)
    const BigInt n_minus_3 = *this - BigInt(3);
    for (size_t i = 0; i < rounds; ++i) {
        // Pick random witness a in [2, n-2]
        BigInt a = BigInt::randomBelow(n_minus_3) + BigInt(2);

        if (!millerRabinWitness(a, d, s, *this)) {
            return false;
        }
    }

    return true;
}

BigInt BigInt::randomBits(size_t numBits)
{
    if (numBits == 0) return BigInt();

    // Use random_device for seeding and mt19937_64 for generation
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint32_t> dist32(0, UINT32_MAX);

    // Calculate number of 32-bit words needed
    size_t numWords = (numBits + 31) / 32;
    std::vector<uint32_t> digits(numWords);

    // Fill with random data
    for (size_t i = 0; i < numWords; ++i) {
        digits[i] = dist32(gen);
    }

    // Mask the top word to get exactly numBits
    size_t topBits = numBits % 32;
    if (topBits == 0) topBits = 32;

    // Set the top bit to ensure we have exactly numBits
    uint32_t mask = (1U << topBits) - 1;
    digits.back() &= mask;
    digits.back() |= (1U << (topBits - 1));  // Set MSB

    // Construct BigInt from little-endian digits
    BigInt result;
    result.digits_ = std::move(digits);
    result.positive_ = true;
    return result;
}

BigInt BigInt::randomBelow(const BigInt& max)
{
    if (max <= BigInt(1)) return BigInt();

    size_t bits = max.bitLength();

    // Rejection sampling to get uniform distribution
    BigInt result;
    do {
        result = randomBits(bits);
    } while (result >= max);

    return result;
}

BigInt BigInt::randomPrime(size_t numBits)
{
    if (numBits < 2) return BigInt(2);
    if (numBits == 2) {
        // Use proper RNG instead of rand()
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<int> dist(0, 1);
        return BigInt(dist(gen) ? 2 : 3);
    }

    while (true) {
        // Generate random odd number with exactly numBits bits
        BigInt candidate = randomBits(numBits);

        // Make it odd
        if (candidate.isEven()) {
            candidate.digits_[0] |= 1;
        }

        // Quick rejection by small primes
        bool divisible = false;
        for (uint32_t p : kSmallPrimes) {
            if (candidate == BigInt(p)) {
                return candidate;  // It's a small prime itself
            }
            if ((candidate % BigInt(p)).isZero()) {
                divisible = true;
                break;
            }
        }

        if (divisible) continue;

        // Full primality test
        if (candidate.isProbablePrime()) {
            return candidate;
        }
    }
}

BigInt BigInt::nextPrime() const
{
    if (*this <= BigInt(2)) return BigInt(2);
    if (*this == BigInt(3)) return BigInt(3);

    // Start with odd number >= *this
    BigInt candidate = *this;
    if (candidate.isEven()) {
        candidate += BigInt(1);
    }

    // By Bertrand's postulate, there's always a prime between n and 2n.
    // We use a generous iteration limit based on the prime gap estimate.
    // For n-bit numbers, expected gap is O(n), so we allow O(n^2) iterations.
    const size_t bits = candidate.bitLength();
    const size_t maxIterations = std::max(static_cast<size_t>(1000000), bits * bits * 100);
    size_t iterations = 0;

    while (iterations < maxIterations) {
        ++iterations;

        // Quick rejection by small primes
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

        candidate += BigInt(2);  // Next odd number
    }

    throw std::runtime_error("nextPrime: exceeded maximum iterations");
}

std::ostream& operator<<(std::ostream& out, const BigInt& value) {
  std::string str = value.toStdString(kDefaultOutputBase);
  out << str;
  return out;
}



std::istream& operator >> (std::istream& in, BigInt& bigInt)
{
    std::string bigNumberString;
    in >> bigNumberString;
    bigInt = BigInt(bigNumberString);
    return in;
}

BigInt BarrettReduction(const BigInt& dividend, const BigInt& divisor, const BigInt& mu)
{
    BigInt remainder = dividend - ((dividend.shiftDigitsToLow(divisor.digits_.size() - 1) * mu).shiftDigitsToLow(divisor.digits_.size() + 1) * divisor);
    while(remainder >= divisor)
    {
        remainder -= divisor;
    }
    return remainder;
}

std::string BigInt::toStdString(const uint8_t base) const
{
    std::stringstream bigNumberStringStream;
    std::string bigNumberString;
    if(digits_.empty() || !(*this))
    {
        return "0";
    }
    if(!positive_)
    {
        bigNumberStringStream << '-';
    }
    if (base == kBaseBinary) {
      for (auto it = digits_.crbegin(); it != digits_.crend(); ++it) {
        bigNumberStringStream << std::bitset<sizeof(uint32_t) * 8>(*it);
      }
    } else if (base == kBaseHexadecimal) {
      for (auto it = digits_.crbegin(); it != digits_.crend(); ++it) {
        bigNumberStringStream << std::hex << std::setw(8) << std::setfill('0')
                              << *it;
      }
    } else {  // base == kBaseDecimal
      const BigInt decimal_repr = toBigIntDec();
      for (auto it = decimal_repr.digits_.crbegin();
           it != decimal_repr.digits_.crend(); ++it) {
        bigNumberStringStream << std::dec << std::setw(9) << std::setfill('0')
                              << *it;
      }
    }
    bigNumberString = bigNumberStringStream.str();
    bigNumberString.erase(positive_ ? 0 : 1, bigNumberString.find_first_not_of("-0") - (positive_ ? 0 : 1));
    return bigNumberString;
}

std::vector<uint32_t> BigInt::toStdVectorUint32_t() const {
  std::vector<uint32_t> result = digits_;
  std::reverse(result.begin(), result.end());
  return result;
}

std::vector<uint8_t> BigInt::toStdVectorUint8_t() const
{
    std::vector<uint8_t> stdVectorUint8_t;
    size_t numberOfBytes = byteLength();
    stdVectorUint8_t.reserve(numberOfBytes);
    std::vector<uint32_t>::const_reverse_iterator iterator = digits_.crbegin();
    if(numberOfBytes % sizeof(uint32_t))
    {
        for(size_t indexFirstBytes = 0; indexFirstBytes < numberOfBytes % sizeof(uint32_t); ++indexFirstBytes)
        {
            stdVectorUint8_t.emplace_back(static_cast<uint8_t>(*iterator >> ((numberOfBytes % sizeof(uint32_t) - indexFirstBytes - 1) * 8)));
        }
        ++iterator;
    }
    while(iterator != digits_.crend())
    {
        for(uint8_t indexByte = 0; indexByte < sizeof(uint32_t); ++indexByte)
        {
            stdVectorUint8_t.emplace_back(static_cast<uint8_t>(*iterator >> ((sizeof(uint32_t) - indexByte - 1) * 8)));
        }
        ++iterator;
    }
    return stdVectorUint8_t;
}

BigInt::operator uint64_t() const noexcept
{
    if(digits_.empty())
    {
        return 0;
    }
    if(digits_.size() >= 2)
    {
        return (static_cast<uint64_t>(*std::next(digits_.cbegin())) << 32) | static_cast<uint64_t>(*digits_.cbegin());
    }
    return digits_.front();
}

BigInt::operator uint32_t() const noexcept
{
    return digits_.empty() ? 0 : digits_.front();
}

BigInt::operator uint16_t() const noexcept
{
    return digits_.empty() ? static_cast<uint16_t>(0) : static_cast<uint16_t>(digits_.front());
}

BigInt::operator uint8_t() const noexcept
{
    return digits_.empty() ? static_cast<uint8_t>(0) : static_cast<uint8_t>(digits_.front());
}

BigInt::operator int64_t() const noexcept
{
    const auto magnitude = static_cast<uint64_t>(*this);
    return positive_ ? static_cast<int64_t>(magnitude) : -static_cast<int64_t>(magnitude);
}

BigInt::operator int32_t() const noexcept
{
    const auto magnitude = static_cast<uint32_t>(*this);
    return positive_ ? static_cast<int32_t>(magnitude) : -static_cast<int32_t>(magnitude);
}

BigInt::operator int16_t() const noexcept
{
    const auto magnitude = static_cast<uint16_t>(*this);
    return positive_ ? static_cast<int16_t>(magnitude) : -static_cast<int16_t>(magnitude);
}

BigInt::operator int8_t() const noexcept
{
    const auto magnitude = static_cast<uint8_t>(*this);
    return positive_ ? static_cast<int8_t>(magnitude) : -static_cast<int8_t>(magnitude);
}

BigInt::operator bool() const noexcept
{
    for(std::vector<uint32_t>::const_iterator iteratordigits_ = digits_.cbegin(); iteratordigits_ != digits_.cend(); ++iteratordigits_)
    {
        if(*iteratordigits_ != 0)
        {
            return true;
        }
    }
    return false;
}

size_t BigInt::bitLength() const noexcept
{
    if(!(*this))
    {
        return 1;
    }
    size_t len = (digits_.size() - 1) * sizeof(uint32_t) * 8;
    uint32_t highOrderDigit = digits_.back();
    uint8_t highOrderBits = 0;
    while(highOrderDigit)
    {
        highOrderDigit >>= 1;
        ++highOrderBits;
    }
    len += highOrderBits;
    return len;
}

size_t BigInt::byteLength() const noexcept
{
    if(!(*this))
    {
        return 1;
    }
    size_t len = (digits_.size() - 1) * sizeof(uint32_t);
    uint32_t highOrderDigit = digits_.back();
    uint8_t highOrderBytes = 0;
    while(highOrderDigit)
    {
        highOrderDigit >>= 8;
        ++highOrderBytes;
    }
    len += highOrderBytes;
    return len;
}

size_t BigInt::digitCount() const noexcept
{
    return digits_.size();
}

bool BigInt::isZero() const noexcept
{
    return digits_.empty() || (digits_.size() == 1 && digits_[0] == 0);
}

bool BigInt::isEven() const noexcept
{
    return digits_.empty() || !(digits_.front() & 1);
}

bool BigInt::isOdd() const noexcept
{
    return !digits_.empty() && (digits_.front() & 1);
}

bool BigInt::isPositive() const noexcept
{
    return positive_ && !isZero();
}

bool BigInt::isNegative() const noexcept
{
    return !positive_ && !isZero();
}

void BigInt::alignTo(BigInt& aligned)
{
    if(digits_.size() > aligned.digits_.size())
    {
        aligned.digits_.reserve(digits_.size());
        aligned.digits_.resize(digits_.size(), 0);
    }
    else if(aligned.digits_.size() > digits_.size())
    {
        digits_.reserve(aligned.digits_.size());
        digits_.resize(aligned.digits_.size(), 0);
    }
}

void BigInt::deleteZeroHighOrderDigit()
{
    while(digits_.size() > 1 && !digits_.back())
    {
        digits_.pop_back();
    }
    // Ensure we always have at least one digit (normalized zero representation)
    if(digits_.empty())
    {
        digits_.push_back(0);
    }
}

BigInt BigInt::shiftDigitsToHigh(const size_t shift) const
{
    BigInt shifted = *this;
    shifted.digits_.insert(shifted.digits_.begin(), shift, 0);
    return shifted;
}

BigInt BigInt::shiftDigitsToLow(const size_t shift) const
{
    BigInt shifted = *this;
    if(shifted.digits_.size() > shift)
    {
        shifted.digits_.erase(shifted.digits_.begin(), shifted.digits_.begin() + static_cast<long>(shift));
    }
    else // shifted.digits_.size() <= shift
    {
        shifted.digits_.shrink_to_fit();
        shifted.digits_.reserve(1);
        shifted.digits_.emplace_back(0);
        shifted.positive_ = true;
    }
    return shifted;
}

BigInt BigInt::multiplySchoolbook(const BigInt& other) const
{
    // Handle zero cases
    if (digits_.empty() || other.digits_.empty()) {
        return BigInt();
    }
    if (isZero() || other.isZero()) {
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
                               static_cast<uint64_t>(result.digits_[i + j]) + carry;
            result.digits_[i + j] = static_cast<uint32_t>(product & UINT32_MAX);
            carry = product >> 32;
        }
        result.digits_[i + n] = static_cast<uint32_t>(carry);
    }

    result.deleteZeroHighOrderDigit();
    result.positive_ = true;  // Sign handled by caller
    return result;
}

BigInt BigInt::multiplyKaratsuba(const BigInt& other) const
{
    // Handle zero cases
    if (isZero() || other.isZero()) {
        return BigInt();
    }

    const size_t m = digits_.size();
    const size_t n = other.digits_.size();

    // Base case: use schoolbook for small numbers
    if (m < kKaratsubaThreshold || n < kKaratsubaThreshold) {
        return multiplySchoolbook(other);
    }

    // Split point: half of the larger operand
    const size_t half = (std::max(m, n) + 1) / 2;

    // Split this = high1 * B^half + low1
    BigInt low1, high1;
    if (m <= half) {
        low1 = *this;
        low1.positive_ = true;
        // high1 is zero (default)
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
        // high2 is zero (default)
    } else {
        low2.digits_.assign(other.digits_.begin(), other.digits_.begin() + static_cast<long>(half));
        low2.positive_ = true;
        low2.deleteZeroHighOrderDigit();
        high2.digits_.assign(other.digits_.begin() + static_cast<long>(half), other.digits_.end());
        high2.positive_ = true;
        high2.deleteZeroHighOrderDigit();
    }

    // Karatsuba: compute 3 products instead of 4
    // z0 = low1 * low2
    // z2 = high1 * high2
    // z1 = (low1 + high1) * (low2 + high2) - z0 - z2
    BigInt z0 = low1.multiplyKaratsuba(low2);
    BigInt z2 = high1.multiplyKaratsuba(high2);

    BigInt sum1 = abs(low1) + abs(high1);
    BigInt sum2 = abs(low2) + abs(high2);
    BigInt z1 = sum1.multiplyKaratsuba(sum2) - z0 - z2;

    // Result = z2 * B^(2*half) + z1 * B^half + z0
    BigInt result = z2.shiftDigitsToHigh(2 * half) + z1.shiftDigitsToHigh(half) + z0;
    result.positive_ = true;  // Sign handled by caller
    return result;
}

BigInt BigInt::toBigIntDec() const
{
    const BigInt kBasisCalcSysDec(1000000000);
    BigInt bigNumber = abs(*this);  // Work with absolute value
    BigInt bigNumberDec;
    std::pair<BigInt, BigInt> BigNumberDivModkBasisCalcSysDec;
    bigNumberDec.positive_ = positive_;
    bigNumberDec.digits_.reserve(digits_.size() + 1);
    while(bigNumber)
    {
        BigNumberDivModkBasisCalcSysDec = bigNumber.DivMod(kBasisCalcSysDec);
        bigNumberDec.digits_.emplace_back(BigNumberDivModkBasisCalcSysDec.second.digits_.front());
        bigNumber = BigNumberDivModkBasisCalcSysDec.first;
    }
    return bigNumberDec;
}

} // namespace bigint



