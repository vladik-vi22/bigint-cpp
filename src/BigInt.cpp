#include <bigint/BigInt.hpp>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace bigint {

namespace {

// Internal constants (not exposed in public API)
constexpr uint64_t kBasisCalcSys = 1ULL << 32;  // 2^32 for carry calculations
constexpr uint32_t kBasisCalcDec = 1000000000;  // 10^9 for decimal conversion
constexpr uint8_t kDecimalCellSize = 9;         // Digits per decimal cell

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

// Define external constants
namespace constants {

const BigInt kZero(static_cast<uint32_t>(0));
const BigInt kOne(static_cast<uint32_t>(1));
const BigInt kTwo(static_cast<uint32_t>(2));
const BigInt kThree(static_cast<uint32_t>(3));
const BigInt kFour(static_cast<uint32_t>(4));
const BigInt kFive(static_cast<uint32_t>(5));
const BigInt kEight(static_cast<uint32_t>(8));

}  // namespace constants

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
  const uint8_t cell_size = base == kBaseHexadecimal
                                ? (sizeof(uint32_t) * 2)
                                : (sizeof(uint32_t) * 8);
  if (!str.empty() && str[0] == '-') {
    positive_ = false;
    str.erase(0, 1);
  } else {
    positive_ = true;
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

BigInt::BigInt(const std::vector<uint32_t>& vec, const bool is_positive_) {
  digits_ = vec;
  std::reverse(digits_.begin(), digits_.end());
  deleteZeroHighOrderDigit();
  positive_ = is_positive_;
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

BigInt::BigInt(const uint64_t value, const bool is_positive_) {
  digits_.emplace_back(static_cast<uint32_t>(value & UINT32_MAX));
  digits_.emplace_back(static_cast<uint32_t>(value >> 32));
  positive_ = is_positive_;
}

BigInt::BigInt(const uint32_t value, const bool is_positive_) {
  digits_.emplace_back(value);
  positive_ = is_positive_;
}

BigInt::BigInt(const int64_t value) {
  const auto abs_value = static_cast<uint64_t>(std::abs(value));
  digits_.emplace_back(static_cast<uint32_t>(abs_value & UINT32_MAX));
  digits_.emplace_back(static_cast<uint32_t>(abs_value >> 32));
  positive_ = (value >= 0);
}

BigInt::BigInt(const int32_t value) {
  digits_.emplace_back(static_cast<uint32_t>(std::abs(value)));
  positive_ = (value >= 0);
}

BigInt& BigInt::operator=(const BigInt& other) {
  digits_ = other.digits_;
  positive_ = other.positive_;
  return *this;
}

BigInt& BigInt::operator=(BigInt&& other) noexcept {
  digits_ = std::move(other.digits_);
  positive_ = other.positive_;
  other.positive_ = true;
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
    *this = *this + augend;
    return *this;
}

BigInt& BigInt::operator ++()
{
    *this += constants::kOne;
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
    *this = *this - subtrahend;
    return *this;
}

BigInt& BigInt::operator -- ()
{
    *this -= constants::kOne;
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
    *this = *this * multiplier;
    return *this;
}

BigInt BigInt::operator * (const BigInt& multiplier) const
{
    BigInt product(0);
    uint32_t shift = 0;
    product.digits_.reserve(digits_.size() + multiplier.digits_.size());
    for(std::vector<uint32_t>::const_iterator iteratorMultiplier = multiplier.digits_.cbegin(); iteratorMultiplier != multiplier.digits_.cend(); ++iteratorMultiplier, ++shift)
    {
        product += (*this * *iteratorMultiplier).shiftDigitsToHigh(shift);
    }
    product.positive_ = positive_ == multiplier.positive_;
    return product;
}

BigInt& BigInt::operator *= (const BigInt& multiplier)
{
    *this = *this * multiplier;
    return *this;
}

std::pair<BigInt, BigInt> BigInt::DivMod(const BigInt& divisor) const
{
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
        fraction += constants::kOne << differenceRemainderNDivisorbitLength; // 1 << n = 2 ^ n
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
        return constants::kZero;
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

size_t log2(const BigInt& antilogarithm)
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
    if(!divisor || !isCoprime(dividend, divisor))
    {
        return constants::kZero;
    }
    BigInt divisor_copy(divisor);
    BigInt fraction;
    BigInt x0(0);
    BigInt x1(1);
    BigInt x_temp;
    while(dividend > constants::kOne)
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

bool congruencemod(const BigInt& dividend1, const BigInt& dividend2, const BigInt divisor)
{
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
    return gcd(bigInt1, bigInt2) == constants::kOne;
}

int8_t symbolJacobi(BigInt bigInt1, BigInt bigInt2)
{
    if(!isCoprime(bigInt1, bigInt2))
    {
        return 0;
    }
    int8_t symbolJacobi = 1;
    size_t iterator = 0;
    BigInt bigInt3;
    if(!bigInt1.positive_)
    {
        bigInt1.positive_ = true;
        if(bigInt2 % constants::kFour == constants::kThree)
        {
            symbolJacobi = -symbolJacobi;
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
            if(bigInt2 % constants::kEight == constants::kThree || bigInt2 % constants::kEight == constants::kFive)
            {
                symbolJacobi = -symbolJacobi;
            }
        }
        if(bigInt1 % constants::kFour == constants::kThree && bigInt2 % constants::kFour == constants::kThree)
        {
            symbolJacobi = -symbolJacobi;
        }
        bigInt3 = bigInt1;
        bigInt1 = bigInt2 % bigInt3;
        bigInt2 = bigInt3;
    }
    return symbolJacobi;
}

BigInt BigInt::operator ~() const
{
    return -*this - constants::kOne;
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
        return -(~*this | ~rightBitwiseAND) - constants::kOne;
    }
    else if(positive_ && !rightBitwiseAND.positive_)
    {
        return (*this | ~rightBitwiseAND) + rightBitwiseAND + constants::kOne;
    }
    else // !positive_ && rightBitwiseAnd.positive_
    {
        return (~*this | rightBitwiseAND) + *this + constants::kOne;
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
        return -(~*this & ~rightBitwiseOR) - constants::kOne;
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
    BigInt shifted;
    shifted.positive_ = positive_;
    if(shift < 32)
    {
        shifted.digits_.reserve(digits_.size() + 1);
        uint32_t carry = 0;
        uint32_t shifted_temp = 0;
        for(std::vector<uint32_t>::const_iterator iteratorShifting = digits_.cbegin(); iteratorShifting != digits_.cend(); ++iteratorShifting)
        {
            shifted_temp = (*iteratorShifting << shift) | carry;
            carry = *iteratorShifting >> (32 - shift);
            shifted.digits_.emplace_back(shifted_temp);
        }
        if(carry)
        {
            shifted.digits_.emplace_back(carry);
        }
        return shifted;
    }
    else // shift >= 32
    {
        shifted = *this;
        for(size_t indexShift = 0; indexShift < shift / (32 - 1); ++indexShift)
        {
            shifted <<= (32 - 1);
        }
        shifted <<= (shift % (32 - 1));
        return shifted;
    }
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
    BigInt shifted;
    shifted.positive_ = positive_;
    if(shift < 32)
    {
        shifted.digits_.reserve(digits_.size());
        uint32_t carry = 0;
        uint32_t shifted_temp = 0;
        const uint32_t mask = static_cast<uint32_t>((1ULL << shift) - 1);  // Mask for lowest 'shift' bits
        // Process from high to low, building result in reverse order
        for(std::vector<uint32_t>::const_reverse_iterator iteratorShifting = digits_.crbegin(); iteratorShifting != digits_.crend(); ++iteratorShifting)
        {
            shifted_temp = (*iteratorShifting >> shift) | carry;
            carry = (*iteratorShifting & mask) << (32 - shift);
            shifted.digits_.emplace_back(shifted_temp);
        }
        // Reverse to get correct order (low to high)
        std::reverse(shifted.digits_.begin(), shifted.digits_.end());
        shifted.deleteZeroHighOrderDigit();
        return shifted;
    }
    else // shift >= 32
    {
        shifted = *this;
        for(uint32_t indexShift = 0; indexShift < shift / (32 - 1); ++indexShift)
        {
            shifted >>= (32 - 1);
        }
        shifted >>= (shift % (32 - 1));
        shifted.deleteZeroHighOrderDigit();
        return shifted;
    }
}

BigInt& BigInt::operator >>= (const size_t shift)
{
    *this = *this >> shift;
    return *this;
}

BigInt BigInt::leftCircularShift(const size_t shift) const
{
    const BigInt mask(--(constants::kOne << bitLength()));
    return (((*this << shift) | (*this >> (bitLength() - shift))) & mask);
}

BigInt BigInt::rightCircularShift(const size_t shift) const
{
    const BigInt mask(--(constants::kOne << bitLength()));
    return (((*this >> shift) | (*this << (bitLength() - shift))) & mask);
}

bool BigInt::operator !() const noexcept
{
    {
        for(std::vector<uint32_t>::const_iterator iteratordigits_ = digits_.cbegin(); iteratordigits_ != digits_.cend(); ++iteratordigits_)
        {
            if(*iteratordigits_ != 0)
            {
                return false;
            }
        }
        return true;
    }
}

bool BigInt::operator && (const BigInt& rightAND) const noexcept
{
    return static_cast<bool>(*this) && static_cast<bool>(rightAND);
}

bool BigInt::operator || (const BigInt& rightOR) const noexcept
{
    return static_cast<bool>(*this) || static_cast<bool>(rightOR);
}

bool BigInt::operator == (const BigInt& rightComparable) const noexcept
{
    // Handle -0 == +0 case
    if (isZero() && rightComparable.isZero()) return true;
    return (digits_ == rightComparable.digits_ && positive_ == rightComparable.positive_);
}

bool BigInt::operator != (const BigInt& rightComparable) const noexcept
{
    // Handle -0 == +0 case
    if (isZero() && rightComparable.isZero()) return false;
    return (positive_ != rightComparable.positive_ || digits_ != rightComparable.digits_);
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

bool BigInt::operator < (const BigInt& rightComparable) const noexcept
{
    // Handle zero cases: 0 is never less than 0
    if (isZero() && rightComparable.isZero()) return false;

    if (positive_ && rightComparable.positive_)
    {
        return compareMagnitude(rightComparable) < 0;
    }
    else if (positive_ && !rightComparable.positive_)
    {
        return false;
    }
    else if (!positive_ && rightComparable.positive_)
    {
        return true;
    }
    else // !positive_ && !rightComparable.positive_
    {
        // For negative numbers: -5 < -3 means |5| > |3|
        return compareMagnitude(rightComparable) > 0;
    }
}

bool BigInt::operator > (const BigInt& rightComparable) const noexcept
{
    // Handle zero cases: 0 is never greater than 0
    if (isZero() && rightComparable.isZero()) return false;

    if (positive_ && rightComparable.positive_)
    {
        return compareMagnitude(rightComparable) > 0;
    }
    else if (positive_ && !rightComparable.positive_)
    {
        return true;
    }
    else if (!positive_ && rightComparable.positive_)
    {
        return false;
    }
    else // !positive_ && !rightComparable.positive_
    {
        // For negative numbers: -3 > -5 means |3| < |5|
        return compareMagnitude(rightComparable) < 0;
    }
}

bool BigInt::operator <= (const BigInt& rightComparable) const noexcept
{
    return !(*this > rightComparable);
}

bool BigInt::operator >= (const BigInt& rightComparable) const noexcept
{
    return !(*this < rightComparable);
}

BigInt abs(const BigInt& bigInt)
{
    BigInt absolute(bigInt);
    absolute.positive_ = true;
    return absolute;
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

bool BigInt::isZero() const noexcept
{
    return digits_.empty() || (digits_.size() == 1 && digits_[0] == 0);
}

bool BigInt::isEven() const noexcept
{
    return !(digits_.front() & 1);
}

bool BigInt::isOdd() const noexcept
{
    return digits_.front() & 1;
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

BigInt BigInt::toBigIntDec() const
{
    const BigInt kBasisCalcSysDec(1000000000);
    BigInt bigNumber(*this);
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



