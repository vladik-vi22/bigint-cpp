#include <bigint/BigInt.hpp>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace bigint {

BigInt::BigInt()
{
}

BigInt::BigInt(const BigInt& bigInt)
{
    vectorUint32_t = bigInt.vectorUint32_t;
    positive = bigInt.positive;
}

BigInt::BigInt(std::string bigIntString, const uint8_t base)
{
    const uint8_t sizeOfCell = base == baseHexadecimal ? (sizeof(uint32_t) * 2) : (sizeof(uint32_t) * 8);
    if(bigIntString[0] == '-')
    {
        positive = false;
        bigIntString.erase(0, 1);
    }
    else
    {
        positive = true;
    }
    if(base == baseDecimal)
    {
        bigIntString = strDec2strBin(bigIntString);
    }
    while(bigIntString.length() % sizeOfCell)
    {
        bigIntString.insert(0, 1, '0');
    }
    size_t sizeOfArr = bigIntString.length() / sizeOfCell;
    vectorUint32_t.reserve(sizeOfArr);
    for(size_t indexVectorUint32_t = 0; indexVectorUint32_t < sizeOfArr; ++indexVectorUint32_t)
    {
        vectorUint32_t.emplace(vectorUint32_t.begin(), static_cast<uint32_t>(std::stoul(bigIntString.substr(indexVectorUint32_t * sizeOfCell, sizeOfCell), nullptr, base == baseHexadecimal ? baseHexadecimal : baseBinary)));
    }
}

BigInt::BigInt(const std::vector<uint32_t>& bigIntVectorUint32_t, const bool isPositive)
{
    vectorUint32_t = bigIntVectorUint32_t;
    std::reverse(vectorUint32_t.begin(), vectorUint32_t.end());
    deleteZeroHighOrderDigit();
    positive = isPositive;
}

BigInt::BigInt(const std::vector<uint16_t>& bigIntVectorUint16_t, const bool isPositive)
{
    vectorUint32_t.reserve(bigIntVectorUint16_t.size() & 1 ? (bigIntVectorUint16_t.size() >> 1) + 1 : bigIntVectorUint16_t.size() >> 1);
    std::vector<uint16_t>::const_reverse_iterator iteratorVectorUint16_t = bigIntVectorUint16_t.crbegin();
    for(size_t indexVectorUint32_t = 0; indexVectorUint32_t < (bigIntVectorUint16_t.size() >> 1); ++indexVectorUint32_t)
    {
        vectorUint32_t.emplace_back(static_cast<uint32_t>(*iteratorVectorUint16_t) |
                                    static_cast<uint32_t>(*(++iteratorVectorUint16_t)) << 16);
        ++iteratorVectorUint16_t;
    }
    if(bigIntVectorUint16_t.size() & 1)
    {
        vectorUint32_t.emplace_back(static_cast<uint32_t>(*iteratorVectorUint16_t));
    }
    deleteZeroHighOrderDigit();
    positive = isPositive;
}

BigInt::BigInt(const std::vector<uint8_t>& bigIntVectorUint8_t, const bool isPositive)
{
    vectorUint32_t.reserve(bigIntVectorUint8_t.size() & 3 ? (bigIntVectorUint8_t.size() >> 2) + 1 : bigIntVectorUint8_t.size() >> 2);
    std::vector<uint8_t>::const_reverse_iterator iteratorVectorUint8_t = bigIntVectorUint8_t.crbegin();
    for(size_t indexVectorUint32_t = 0; indexVectorUint32_t < (bigIntVectorUint8_t.size() >> 2); ++indexVectorUint32_t)
    {
        vectorUint32_t.emplace_back(static_cast<uint32_t>(*iteratorVectorUint8_t) |
                                    static_cast<uint32_t>(*(++iteratorVectorUint8_t)) << 8 |
                                    static_cast<uint32_t>(*(++iteratorVectorUint8_t)) << 16 |
                                    static_cast<uint32_t>(*(++iteratorVectorUint8_t)) << 24);
        ++iteratorVectorUint8_t;
    }
    if((bigIntVectorUint8_t.size() & 3) == 3)
    {
        vectorUint32_t.emplace_back(static_cast<uint32_t>(*iteratorVectorUint8_t) |
                                    static_cast<uint32_t>(*(++iteratorVectorUint8_t)) << 8 |
                                    static_cast<uint32_t>(*(++iteratorVectorUint8_t)) << 16);
    }
    else if((bigIntVectorUint8_t.size() & 3) == 2)
    {
        vectorUint32_t.emplace_back(static_cast<uint32_t>(*iteratorVectorUint8_t) |
                                    static_cast<uint32_t>(*(++iteratorVectorUint8_t)) << 8);
    }
    else if((bigIntVectorUint8_t.size() & 3) == 1)
    {
        vectorUint32_t.emplace_back(static_cast<uint32_t>(*iteratorVectorUint8_t));
    }
    deleteZeroHighOrderDigit();
    positive = isPositive;
}

BigInt::BigInt(const std::vector<bool>& bigIntVectorBool, const bool isPositive)
{
    vectorUint32_t.reserve(bigIntVectorBool.size() & 31 ? (bigIntVectorBool.size() >> 5) + 1 : bigIntVectorBool.size() >> 5);
    uint32_t vectorUint32_t_element;
    std::vector<bool>::const_reverse_iterator iteratorVectorBool = bigIntVectorBool.crbegin();
    for(size_t indexVectorUint32_t = 0; indexVectorUint32_t < (bigIntVectorBool.size() >> 5); ++indexVectorUint32_t)
    {
        vectorUint32_t_element = 0;
        for(uint8_t indexBit = 0; indexBit < 32; ++indexBit)
        {
            vectorUint32_t_element |= static_cast<uint32_t>(*iteratorVectorBool) << indexBit;
            ++iteratorVectorBool;
        }
        vectorUint32_t.emplace_back(vectorUint32_t_element);
    }
    if(bigIntVectorBool.size() & 31)
    {
        vectorUint32_t_element = 0;
        for(uint8_t indexBit = 0; indexBit < (bigIntVectorBool.size() & 31); ++indexBit)
        {
            vectorUint32_t_element |= static_cast<uint32_t>(*iteratorVectorBool) << indexBit;
            ++iteratorVectorBool;
        }
        vectorUint32_t.emplace_back(vectorUint32_t_element);
    }
    deleteZeroHighOrderDigit();
    positive = isPositive;
}

BigInt::BigInt(const uint64_t bigIntUint64_t, const bool isPositive)
{
    vectorUint32_t.emplace_back(bigIntUint64_t & UINT32_MAX);
    vectorUint32_t.emplace_back(bigIntUint64_t >> 32);
    positive = isPositive;
}

BigInt::BigInt(const uint32_t bigIntUint32_t, const bool isPositive)
{
    vectorUint32_t.emplace_back(bigIntUint32_t);
    positive = isPositive;
}

BigInt::BigInt(const int64_t bigIntInt64_t)
{
    vectorUint32_t.emplace_back(std::abs(bigIntInt64_t) & UINT32_MAX);
    vectorUint32_t.emplace_back(static_cast<uint32_t>(std::abs(bigIntInt64_t) >> 32));
    positive = (bigIntInt64_t >= 0);
}

BigInt::BigInt(const int32_t bigIntInt32_t)
{
    vectorUint32_t.emplace_back(static_cast<uint32_t>(std::abs(bigIntInt32_t)));
    positive = (bigIntInt32_t >= 0);
}

BigInt& BigInt::operator = (const BigInt& equal)
{
    vectorUint32_t = equal.vectorUint32_t;
    positive = equal.positive;
    return *this;
}

BigInt BigInt::operator +() const
{
    return *this;
}

BigInt BigInt::operator + (const BigInt& addend) const
{
    if(positive && addend.positive)
    {
        BigInt sum;
        uint32_t carry = 0;
        uint64_t sum_temp;
        bool augendGreater = (vectorUint32_t.size() >= addend.vectorUint32_t.size());
        sum.vectorUint32_t.reserve(augendGreater ? vectorUint32_t.size() + 1 : addend.vectorUint32_t.size() + 1);
        std::vector<uint32_t>::const_iterator iteratorAugend = augendGreater ? vectorUint32_t.cbegin() : addend.vectorUint32_t.cbegin();
        std::vector<uint32_t>::const_iterator iteratorAddend = augendGreater ? addend.vectorUint32_t.cbegin() : vectorUint32_t.cbegin();
        std::vector<uint32_t>::const_iterator iteratorAugendEnd = augendGreater ? vectorUint32_t.cend() : addend.vectorUint32_t.cend();
        std::vector<uint32_t>::const_iterator iteratorAddendEnd = augendGreater ? addend.vectorUint32_t.cend() : vectorUint32_t.cend();
        while(iteratorAddend != iteratorAddendEnd)
        {
            sum_temp = static_cast<uint64_t>(*iteratorAugend) + static_cast<uint64_t>(*iteratorAddend) + static_cast<uint64_t>(carry);
            sum.vectorUint32_t.emplace_back(sum_temp & UINT32_MAX);
            carry = sum_temp >> 32;
            ++iteratorAugend;
            ++iteratorAddend;
        }
        while(iteratorAugend != iteratorAugendEnd)
        {
            sum_temp = static_cast<uint64_t>(*iteratorAugend) + static_cast<uint64_t>(carry);
            sum.vectorUint32_t.emplace_back(sum_temp & UINT32_MAX);
            carry = sum_temp >> 32;
            ++iteratorAugend;
        }
        if(carry)
        {
            sum.vectorUint32_t.emplace_back(carry);
        }
        sum.positive = true;
        return sum;
    }
    else if(positive && !addend.positive)
    {
        return *this - abs(addend);
    }
    else if(!positive && addend.positive)
    {
        return addend - abs(*this);
    }
    else // !positive && !addend.positive
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
    *this += constants::One;
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
    negative.positive = !positive;
    return negative;
}

BigInt BigInt::operator - (const BigInt& subtrahend) const
{
    if(positive && subtrahend.positive)
    {
        if(*this >= subtrahend)
        {
            BigInt difference;
            uint8_t borrow = 0;
            int64_t difference_temp;
            difference.vectorUint32_t.reserve(vectorUint32_t.size());
            std::vector<uint32_t>::const_iterator iteratorMinuend = vectorUint32_t.cbegin();
            std::vector<uint32_t>::const_iterator iteratorSubtrahend = subtrahend.vectorUint32_t.cbegin();
            while(iteratorSubtrahend != subtrahend.vectorUint32_t.cend())
            {
                difference_temp = static_cast<int64_t>(*iteratorMinuend) - static_cast<int64_t>(*iteratorSubtrahend) - static_cast<int64_t>(borrow);
                if(difference_temp >= 0)
                {
                    difference.vectorUint32_t.emplace_back(static_cast<uint32_t>(difference_temp));
                    borrow = 0;
                }
                else // difference_temp < 0
                {
                    difference.vectorUint32_t.emplace_back(static_cast<uint32_t>(difference_temp + static_cast<int64_t>(basisCalcSys)));
                    borrow = 1;
                }
                ++iteratorMinuend;
                ++iteratorSubtrahend;
            }
            while(iteratorMinuend != vectorUint32_t.cend())
            {
                difference_temp = static_cast<int64_t>(*iteratorMinuend) - static_cast<int64_t>(borrow);
                if(difference_temp >= 0)
                {
                    difference.vectorUint32_t.emplace_back(static_cast<uint32_t>(difference_temp));
                    borrow = 0;
                }
                else // difference_temp < 0
                {
                    difference.vectorUint32_t.emplace_back(static_cast<uint32_t>(difference_temp + static_cast<int64_t>(basisCalcSys)));
                    borrow = 1;
                }
                ++iteratorMinuend;
            }
            difference.deleteZeroHighOrderDigit();
            difference.positive = true;
            return difference;
        }
        else // minuend < subtrahend
        {
            return -(subtrahend - *this);
        }
    }
    else if(!positive && subtrahend.positive)
    {
        return -(abs(*this) + subtrahend);
    }
    else if(positive && !subtrahend.positive)
    {
        return *this + abs(subtrahend);
    }
    else // !positive && !subtrahend.positive
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
    *this -= constants::One;
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
    product.vectorUint32_t.reserve(vectorUint32_t.size() + 1);
    for(std::vector<uint32_t>::const_iterator iteratorMultiplicand = vectorUint32_t.cbegin(); iteratorMultiplicand != vectorUint32_t.cend(); ++iteratorMultiplicand)
    {
        product_temp = static_cast<uint64_t>(*iteratorMultiplicand) * static_cast<uint64_t>(multiplier) + static_cast<uint64_t>(carry);
        product.vectorUint32_t.emplace_back(product_temp & UINT32_MAX);
        carry = product_temp >> 32;
    }
    if(carry)
    {
        product.vectorUint32_t.emplace_back(carry);
    }
    product.positive = positive;
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
    product.vectorUint32_t.reserve(vectorUint32_t.size() + multiplier.vectorUint32_t.size());
    for(std::vector<uint32_t>::const_iterator iteratorMultiplier = multiplier.vectorUint32_t.cbegin(); iteratorMultiplier != multiplier.vectorUint32_t.cend(); ++iteratorMultiplier, ++shift)
    {
        product += (*this * *iteratorMultiplier).shiftDigitsToHigh(shift);
    }
    product.positive = positive == multiplier.positive;
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
    fraction.vectorUint32_t.reserve(vectorUint32_t.size());
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
        fraction += constants::One << differenceRemainderNDivisorbitLength; // 1 << n = 2 ^ n
    }
    fraction.positive = positive == divisor.positive;
    remainder.positive = positive;
    if(divisor.positive)
    {
        while(!remainder.positive)
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
    if(!exponent.positive)
    {
        return constants::Zero;
    }
    BigInt power(1);
    power.vectorUint32_t.reserve(base.vectorUint32_t.size() * static_cast<size_t>(exponent));
    for(size_t indexBitExponent = exponent.bitLength() - 1; indexBitExponent > 0; --indexBitExponent)
    {
        if(exponent.vectorUint32_t[indexBitExponent >> 5] & (1 << (indexBitExponent & 31)))
        {
            power *= base;
        }
        power *= power;
    }
    if(exponent.isOdd())
    {
        power *= base;
    }
    if(!base.positive)
    {
        power.positive = exponent.isEven();
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
    power.vectorUint32_t.reserve(divisor.vectorUint32_t.size());
    const BigInt mu = power.shiftDigitsToHigh(divisor.vectorUint32_t.size() * 2) / divisor;
    const uint32_t bitLengthExponent = exponent.bitLength();
    for(size_t indexBitExponent = 0; indexBitExponent < bitLengthExponent; ++indexBitExponent)
    {
        if(exponent.vectorUint32_t[indexBitExponent >> 5] & (1 << (indexBitExponent & 31)))
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
    power.vectorUint32_t.reserve(divisor.vectorUint32_t.size());
    const size_t bitLen = exponent.bitLength();
    for(size_t i = 0; i < bitLen; ++i)
    {
        size_t indexBitExponent = bitLen - 1 - i;
        power = (power * power) % divisor;
        if(exponent.vectorUint32_t[indexBitExponent >> 5] & (1 << (indexBitExponent & 31)))
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
        return constants::Zero;
    }
    BigInt divisor_copy(divisor);
    BigInt fraction;
    BigInt x0(0);
    BigInt x1(1);
    BigInt x_temp;
    while(dividend > constants::One)
    {
        fraction = dividend / divisor_copy;
        x_temp = divisor_copy;
        divisor_copy = dividend % divisor_copy;
        dividend = x_temp;
        x_temp = x0;
        x0 = x1 - (fraction * x0);
        x1 = x_temp;
    }
    if(!x1.positive)
    {
        x1 += divisor;
    }
    return x1;
}

bool congruencemod(const BigInt& dividend1, const BigInt& dividend2, const BigInt divisor)
{
    BigInt remainder1(dividend1 % divisor);
    BigInt remainder2(dividend2 % divisor);
    while(!remainder1.positive)
    {
        remainder1 += divisor;
    }
    while(remainder1 > divisor)
    {
        remainder1 -= divisor;
    }
    while(!remainder2.positive)
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
    return gcd(bigInt1, bigInt2) == constants::One;
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
    if(!bigInt1.positive)
    {
        bigInt1.positive = true;
        if(bigInt2 % constants::Four == constants::Three)
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
            if(bigInt2 % constants::Eight == constants::Three || bigInt2 % constants::Eight == constants::Five)
            {
                symbolJacobi = -symbolJacobi;
            }
        }
        if(bigInt1 % constants::Four == constants::Three && bigInt2 % constants::Four == constants::Three)
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
    return -*this - constants::One;
}

BigInt BigInt::operator & (const BigInt& rightBitwiseAND) const
{
    if(positive && rightBitwiseAND.positive)
    {
        BigInt bitwiseAND;
        bitwiseAND.vectorUint32_t.reserve(std::min(vectorUint32_t.size(), rightBitwiseAND.vectorUint32_t.size()));
        std::vector<uint32_t>::const_iterator iteratorLeftBitwiseAND = vectorUint32_t.cbegin();
        std::vector<uint32_t>::const_iterator iteratorRightBitwiseAND = rightBitwiseAND.vectorUint32_t.cbegin();
        while(iteratorLeftBitwiseAND != vectorUint32_t.cend() && iteratorRightBitwiseAND != rightBitwiseAND.vectorUint32_t.cend())
        {
            bitwiseAND.vectorUint32_t.emplace_back(*iteratorLeftBitwiseAND & *iteratorRightBitwiseAND);
            ++iteratorLeftBitwiseAND;
            ++iteratorRightBitwiseAND;
        }
        bitwiseAND.positive = true;
        return bitwiseAND;
    }
    else if(!positive && !rightBitwiseAND.positive)
    {
        return -(~*this | ~rightBitwiseAND) - constants::One;
    }
    else if(positive && !rightBitwiseAND.positive)
    {
        return (*this | ~rightBitwiseAND) + rightBitwiseAND + constants::One;
    }
    else // !positive && rightBitwiseAnd.positive
    {
        return (~*this | rightBitwiseAND) + *this + constants::One;
    }
}

BigInt& BigInt::operator &= (const BigInt& rightBitwiseAND)
{
    *this = *this & rightBitwiseAND;
    return *this;
}

BigInt BigInt::operator | (const BigInt& rightBitwiseOR) const
{
    if(positive && rightBitwiseOR.positive)
    {
        BigInt bitwiseOR;
        const bool leftBitwiseORGreater = (vectorUint32_t.size() >= rightBitwiseOR.vectorUint32_t.size());
        bitwiseOR.vectorUint32_t.reserve(leftBitwiseORGreater ? vectorUint32_t.size() + 1 : rightBitwiseOR.vectorUint32_t.size() + 1);
        std::vector<uint32_t>::const_iterator iteratorLeftBitwiseOR = leftBitwiseORGreater ? vectorUint32_t.cbegin() : rightBitwiseOR.vectorUint32_t.cbegin();
        std::vector<uint32_t>::const_iterator iteratorRightBitwiseOR = leftBitwiseORGreater ? rightBitwiseOR.vectorUint32_t.cbegin() : vectorUint32_t.cbegin();
        std::vector<uint32_t>::const_iterator iteratorLeftBitwiseOREnd = leftBitwiseORGreater ? vectorUint32_t.cend() : rightBitwiseOR.vectorUint32_t.cend();
        std::vector<uint32_t>::const_iterator iteratorRightBitwiseOREnd = leftBitwiseORGreater ? rightBitwiseOR.vectorUint32_t.cend() : vectorUint32_t.cend();
        while(iteratorRightBitwiseOR != iteratorRightBitwiseOREnd)
        {
            bitwiseOR.vectorUint32_t.emplace_back(*iteratorLeftBitwiseOR | *iteratorRightBitwiseOR);
            ++iteratorLeftBitwiseOR;
            ++iteratorRightBitwiseOR;
        }
        while(iteratorLeftBitwiseOR != iteratorLeftBitwiseOREnd)
        {
            bitwiseOR.vectorUint32_t.emplace_back(*iteratorLeftBitwiseOR);
            ++iteratorLeftBitwiseOR;
        }
        bitwiseOR.positive = true;
        return bitwiseOR;
    }
    else if(!positive && !rightBitwiseOR.positive)
    {
        return -(~*this & ~rightBitwiseOR) - constants::One;
    }
    else if(positive && !rightBitwiseOR.positive)
    {
        return (*this & ~rightBitwiseOR) + rightBitwiseOR;
    }
    else // !positive && rightBitwiseOR.positive
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
    if(positive && rightBitwiseXOR.positive)
    {
        BigInt bitwiseXOR;
        const bool leftBitwiseXORGreater = (vectorUint32_t.size() >= rightBitwiseXOR.vectorUint32_t.size());
        bitwiseXOR.vectorUint32_t.reserve(leftBitwiseXORGreater ? vectorUint32_t.size() + 1 : rightBitwiseXOR.vectorUint32_t.size() + 1);
        std::vector<uint32_t>::const_iterator iteratorLeftBitwiseXOR = leftBitwiseXORGreater ? vectorUint32_t.cbegin() : rightBitwiseXOR.vectorUint32_t.cbegin();
        std::vector<uint32_t>::const_iterator iteratorRightBitwiseXOR = leftBitwiseXORGreater ? rightBitwiseXOR.vectorUint32_t.cbegin() : vectorUint32_t.cbegin();
        std::vector<uint32_t>::const_iterator iteratorLeftBitwiseXOREnd = leftBitwiseXORGreater ? vectorUint32_t.cend() : rightBitwiseXOR.vectorUint32_t.cend();
        std::vector<uint32_t>::const_iterator iteratorRightBitwiseXOREnd = leftBitwiseXORGreater ? rightBitwiseXOR.vectorUint32_t.cend() : vectorUint32_t.cend();
        while(iteratorRightBitwiseXOR != iteratorRightBitwiseXOREnd)
        {
            bitwiseXOR.vectorUint32_t.emplace_back(*iteratorLeftBitwiseXOR ^ *iteratorRightBitwiseXOR);
            ++iteratorLeftBitwiseXOR;
            ++iteratorRightBitwiseXOR;
        }
        while(iteratorLeftBitwiseXOR != iteratorLeftBitwiseXOREnd)
        {
            bitwiseXOR.vectorUint32_t.emplace_back(*iteratorLeftBitwiseXOR);
            ++iteratorLeftBitwiseXOR;
        }
        bitwiseXOR.positive = true;
        return bitwiseXOR;
    }
    else // !positive || !rightBitwiseXOR.positive
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
    if(!shift)
    {
        return *this;
    }
    BigInt shifted;
    shifted.positive = positive;
    if(shift < 32)
    {
        shifted.vectorUint32_t.reserve(vectorUint32_t.size() + 1);
        uint32_t carry = 0;
        uint32_t shifted_temp = 0;
        for(std::vector<uint32_t>::const_iterator iteratorShifting = vectorUint32_t.cbegin(); iteratorShifting != vectorUint32_t.cend(); ++iteratorShifting)
        {
            shifted_temp = (*iteratorShifting << shift) | carry;
            carry = *iteratorShifting >> (32 - shift);
            shifted.vectorUint32_t.emplace_back(shifted_temp);
        }
        if(carry)
        {
            shifted.vectorUint32_t.emplace_back(carry);
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
    if(!shift)
    {
        return *this;
    }
    BigInt shifted;
    shifted.positive = positive;
    if(shift < 32)
    {
        shifted.vectorUint32_t.reserve(vectorUint32_t.size());
        uint32_t carry = 0;
        uint32_t shifted_temp = 0;
        const uint32_t mask = static_cast<uint32_t>((1ULL << shift) - 1);  // Mask for lowest 'shift' bits
        // Process from high to low, building result in reverse order
        for(std::vector<uint32_t>::const_reverse_iterator iteratorShifting = vectorUint32_t.crbegin(); iteratorShifting != vectorUint32_t.crend(); ++iteratorShifting)
        {
            shifted_temp = (*iteratorShifting >> shift) | carry;
            carry = (*iteratorShifting & mask) << (32 - shift);
            shifted.vectorUint32_t.emplace_back(shifted_temp);
        }
        // Reverse to get correct order (low to high)
        std::reverse(shifted.vectorUint32_t.begin(), shifted.vectorUint32_t.end());
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
    const BigInt mask(--(constants::One << bitLength()));
    return (((*this << shift) | (*this >> (bitLength() - shift))) & mask);
}

BigInt BigInt::rightCircularShift(const size_t shift) const
{
    const BigInt mask(--(constants::One << bitLength()));
    return (((*this >> shift) | (*this << (bitLength() - shift))) & mask);
}

bool BigInt::operator !() const
{
    {
        for(std::vector<uint32_t>::const_iterator iteratorVectorUint32_t = vectorUint32_t.cbegin(); iteratorVectorUint32_t != vectorUint32_t.cend(); ++iteratorVectorUint32_t)
        {
            if(*iteratorVectorUint32_t != 0)
            {
                return false;
            }
        }
        return true;
    }
}

bool BigInt::operator && (const BigInt& rightAND) const
{
    return static_cast<bool>(*this) && static_cast<bool>(rightAND);
}

bool BigInt::operator || (const BigInt& rightOR) const
{
    return static_cast<bool>(*this) || static_cast<bool>(rightOR);
}

bool BigInt::operator == (const BigInt& rightComparable) const
{
    return (vectorUint32_t == rightComparable.vectorUint32_t && positive == rightComparable.positive);
}

bool BigInt::operator != (const BigInt& rightComparable) const
{
    return (positive != rightComparable.positive || vectorUint32_t != rightComparable.vectorUint32_t);
}

bool BigInt::operator < (const BigInt& rightComparable) const
{
    if(positive && rightComparable.positive)
    {
        if(vectorUint32_t.size() > rightComparable.vectorUint32_t.size())
        {
            return false;
        }
        else if(vectorUint32_t.size() < rightComparable.vectorUint32_t.size())
        {
            return true;
        }
        else // vectorUint32_t.size() == rightComparable.vectorUint32_t.size()
        {
            for(std::vector<uint32_t>::const_reverse_iterator iteratorLeftComparable = vectorUint32_t.crbegin(), iteratorRightComparable = rightComparable.vectorUint32_t.crbegin() ; iteratorLeftComparable != vectorUint32_t.crend(); ++iteratorLeftComparable, ++iteratorRightComparable)
            {
                if(*iteratorLeftComparable != *iteratorRightComparable)
                {
                    return *iteratorLeftComparable < *iteratorRightComparable;
                }
            }
            return false;
        }
    }
    else if(positive && !rightComparable.positive)
    {
        return false;
    }
    else if(!positive && rightComparable.positive)
    {
        return true;
    }
    else // !positive && !rightComparable.positive
    {
        return abs(*this) > abs(rightComparable);
    }
}

bool BigInt::operator > (const BigInt& rightComparable) const
{
    if(positive && rightComparable.positive)
    {
        if(vectorUint32_t.size() > rightComparable.vectorUint32_t.size())
        {
            return true;
        }
        else if(vectorUint32_t.size() < rightComparable.vectorUint32_t.size())
        {
            return false;
        }
        else // vectorUint32_t.size == rightComparable.vectorUint32_t.size()
        {
            for(std::vector<uint32_t>::const_reverse_iterator iteratorLeftComparable = vectorUint32_t.crbegin(), iteratorRightComparable = rightComparable.vectorUint32_t.crbegin(); iteratorLeftComparable != vectorUint32_t.crend(); ++iteratorLeftComparable, ++iteratorRightComparable)
            {
                if(*iteratorLeftComparable != *iteratorRightComparable)
                {
                    return *iteratorLeftComparable > *iteratorRightComparable;
                }
            }
            return false;
        }
    }
    else if(positive && !rightComparable.positive)
    {
        return true;
    }
    else if(!positive && rightComparable.positive)
    {
        return false;
    }
    else // !positive && !rightComparable.positive
    {
        return abs(*this) < abs(rightComparable);
    }
}

bool BigInt::operator <= (const BigInt& rightComparable) const
{
    return (*this == rightComparable || *this < rightComparable);
}

bool BigInt::operator >= (const BigInt& rightComparable) const
{
    return (*this == rightComparable || *this > rightComparable);
}

BigInt abs(const BigInt& bigInt)
{
    BigInt absolute(bigInt);
    absolute.positive = true;
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
    bigInt1.positive = true;
    bigInt2.positive = true;
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

const BigInt& max(const BigInt& bigInt1, const BigInt& bigInt2)
{
    return bigInt1 > bigInt2 ? bigInt1 : bigInt2;
}

const BigInt& min(const BigInt& bigInt1, const BigInt& bigInt2)
{
    return bigInt1 < bigInt2 ? bigInt1 : bigInt2;
}

std::ostream& operator << (std::ostream& out, const BigInt& bigInt)
{
    std::string bigNumberString = bigInt.toStdString(BigInt::baseOutput);
    out << bigNumberString;
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
    BigInt remainder = dividend - ((dividend.shiftDigitsToLow(divisor.vectorUint32_t.size() - 1) * mu).shiftDigitsToLow(divisor.vectorUint32_t.size() + 1) * divisor);
    while(remainder >= divisor)
    {
        remainder -= divisor;
    }
    return remainder;
}

std::string BigInt::toStdString(const int base) const
{
    std::stringstream bigNumberStringStream;
    std::string bigNumberString;
    if(vectorUint32_t.empty() || !(*this))
    {
        return "0";
    }
    if(!positive)
    {
        bigNumberStringStream << '-';
    }
    if(base == baseBinary)
    {
        for(std::vector<uint32_t>::const_reverse_iterator iterator = vectorUint32_t.crbegin(); iterator != vectorUint32_t.crend(); ++iterator)
        {
            bigNumberStringStream << std::bitset<sizeof(uint32_t) * 8>(*iterator);
        }
    }
    else if(base == baseHexadecimal)
    {
        for(std::vector<uint32_t>::const_reverse_iterator iterator = vectorUint32_t.crbegin(); iterator != vectorUint32_t.crend(); ++iterator)
        {
            bigNumberStringStream << std::hex << std::setw(8) << std::setfill('0') << *iterator;
        }
    }
    else // base == baseDecimal
    {
        const BigInt bigNumberDec = toBigIntDec();
        for(std::vector<uint32_t>::const_reverse_iterator iterator = bigNumberDec.vectorUint32_t.crbegin(); iterator != bigNumberDec.vectorUint32_t.crend(); ++iterator)
        {
            bigNumberStringStream << std::dec << std::setw(9) << std::setfill('0') << *iterator;
        }
    }
    bigNumberString = bigNumberStringStream.str();
    bigNumberString.erase(positive ? 0 : 1, bigNumberString.find_first_not_of("-0") - (positive ? 0 : 1));
    return bigNumberString;
}

std::vector<uint32_t> BigInt::toStdVectorUint32_t() const
{
    std::vector<uint32_t> stdVectorUint32_t = vectorUint32_t;
    std::reverse(stdVectorUint32_t.begin(), stdVectorUint32_t.end());
    return stdVectorUint32_t;
}

std::vector<uint8_t> BigInt::toStdVectorUint8_t() const
{
    std::vector<uint8_t> stdVectorUint8_t;
    size_t numberOfBytes = byteLength();
    stdVectorUint8_t.reserve(numberOfBytes);
    std::vector<uint32_t>::const_reverse_iterator iterator = vectorUint32_t.crbegin();
    if(numberOfBytes % sizeof(uint32_t))
    {
        for(size_t indexFirstBytes = 0; indexFirstBytes < numberOfBytes % sizeof(uint32_t); ++indexFirstBytes)
        {
            stdVectorUint8_t.emplace_back(static_cast<uint8_t>(*iterator >> ((numberOfBytes % sizeof(uint32_t) - indexFirstBytes - 1) * 8)));
        }
        ++iterator;
    }
    while(iterator != vectorUint32_t.crend())
    {
        for(uint8_t indexByte = 0; indexByte < sizeof(uint32_t); ++indexByte)
        {
            stdVectorUint8_t.emplace_back(static_cast<uint8_t>(*iterator >> ((sizeof(uint32_t) - indexByte - 1) * 8)));
        }
        ++iterator;
    }
    return stdVectorUint8_t;
}

BigInt::operator uint64_t() const
{
    if(vectorUint32_t.size() >= 2)
    {
        return (static_cast<uint64_t>(*std::next(vectorUint32_t.cbegin())) << 32) | static_cast<uint64_t>(*vectorUint32_t.cbegin());
    }
    else
    {
        return vectorUint32_t.front();
    }
}

BigInt::operator uint32_t() const
{
    return vectorUint32_t.front();
}

BigInt::operator uint16_t() const
{
    return static_cast<uint16_t>(vectorUint32_t.front());
}

BigInt::operator uint8_t() const
{
    return static_cast<uint8_t>(vectorUint32_t.front());
}

BigInt::operator bool() const
{
    for(std::vector<uint32_t>::const_iterator iteratorVectorUint32_t = vectorUint32_t.cbegin(); iteratorVectorUint32_t != vectorUint32_t.cend(); ++iteratorVectorUint32_t)
    {
        if(*iteratorVectorUint32_t != 0)
        {
            return true;
        }
    }
    return false;
}

size_t BigInt::bitLength() const
{
    if(!(*this))
    {
        return 1;
    }
    size_t len = (vectorUint32_t.size() - 1) * sizeof(uint32_t) * 8;
    uint32_t highOrderDigit = vectorUint32_t.back();
    uint8_t highOrderBits = 0;
    while(highOrderDigit)
    {
        highOrderDigit >>= 1;
        ++highOrderBits;
    }
    len += highOrderBits;
    return len;
}

size_t BigInt::byteLength() const
{
    if(!(*this))
    {
        return 1;
    }
    size_t len = (vectorUint32_t.size() - 1) * sizeof(uint32_t);
    uint32_t highOrderDigit = vectorUint32_t.back();
    uint8_t highOrderBytes = 0;
    while(highOrderDigit)
    {
        highOrderDigit >>= 8;
        ++highOrderBytes;
    }
    len += highOrderBytes;
    return len;
}

bool BigInt::isZero() const
{
    return vectorUint32_t.empty() || (vectorUint32_t.size() == 1 && vectorUint32_t[0] == 0);
}

bool BigInt::isEven() const
{
    return !(vectorUint32_t.front() & 1);
}

bool BigInt::isOdd() const
{
    return vectorUint32_t.front() & 1;
}

bool BigInt::isPositive() const
{
    return positive;
}

bool BigInt::isNegative() const
{
    return !positive;
}

void BigInt::alignTo(BigInt& aligned)
{
    if(vectorUint32_t.size() > aligned.vectorUint32_t.size())
    {
        aligned.vectorUint32_t.reserve(vectorUint32_t.size());
        aligned.vectorUint32_t.resize(vectorUint32_t.size(), 0);
    }
    else if(aligned.vectorUint32_t.size() > vectorUint32_t.size())
    {
        vectorUint32_t.reserve(aligned.vectorUint32_t.size());
        vectorUint32_t.resize(aligned.vectorUint32_t.size(), 0);
    }
}

void BigInt::deleteZeroHighOrderDigit()
{
    while(!vectorUint32_t.back() && vectorUint32_t.size() > 1)
    {
        vectorUint32_t.pop_back();
    }
}

BigInt BigInt::shiftDigitsToHigh(const size_t shift) const
{
    BigInt shifted = *this;
    shifted.vectorUint32_t.insert(shifted.vectorUint32_t.begin(), shift, 0);
    return shifted;
}

BigInt BigInt::shiftDigitsToLow(const size_t shift) const
{
    BigInt shifted = *this;
    if(shifted.vectorUint32_t.size() > shift)
    {
        shifted.vectorUint32_t.erase(shifted.vectorUint32_t.begin(), shifted.vectorUint32_t.begin() + static_cast<long>(shift));
    }
    else // shifted.vectorUint32_t.size() <= shift
    {
        shifted.vectorUint32_t.shrink_to_fit();
        shifted.vectorUint32_t.reserve(1);
        shifted.vectorUint32_t.emplace_back(0);
        shifted.positive = true;
    }
    return shifted;
}

BigInt BigInt::toBigIntDec() const
{
    const BigInt basisCalcSysDec(1000000000);
    BigInt bigNumber(*this);
    BigInt bigNumberDec;
    std::pair<BigInt, BigInt> BigNumberDivModBasisCalcSysDec;
    bigNumberDec.positive = positive;
    bigNumberDec.vectorUint32_t.reserve(vectorUint32_t.size() + 1);
    while(bigNumber)
    {
        BigNumberDivModBasisCalcSysDec = bigNumber.DivMod(basisCalcSysDec);
        bigNumberDec.vectorUint32_t.emplace_back(BigNumberDivModBasisCalcSysDec.second.vectorUint32_t.front());
        bigNumber = BigNumberDivModBasisCalcSysDec.first;
    }
    return bigNumberDec;
}

std::string strDec2strBin(std::string strDec)
{
    if(strDec == "0")
    {
        return "0";
    }
    const uint8_t sizeOfCell = 9;
    const uint32_t basisCalc = 1000000000;
    std::string strBin;
    std::vector<uint32_t> vectorUint32_t;
    std::vector<uint32_t> zeroArr;
    uint32_t carryNext;
    uint32_t carryCurrent;
    char charBin;
    while(strDec.length() % sizeOfCell)
    {
        strDec.insert(0, 1, '0');
    }
    size_t sizeOfVector = strDec.length() / sizeOfCell;
    vectorUint32_t.reserve(sizeOfVector);
    for(size_t indexVectorUint32_t = 0; indexVectorUint32_t < sizeOfVector; ++indexVectorUint32_t)
    {
        vectorUint32_t.emplace_back(static_cast<uint32_t>(std::stoul(strDec.substr(indexVectorUint32_t * sizeOfCell, sizeOfCell), nullptr, 10)));
    }
    zeroArr.resize(vectorUint32_t.size(), 0);
    while(vectorUint32_t != zeroArr)
    {
        carryNext = 0;
        for(std::vector<uint32_t>::iterator iteratorShifting = vectorUint32_t.begin(); iteratorShifting != vectorUint32_t.end(); ++iteratorShifting)
        {
            carryCurrent = carryNext;
            carryNext = (*iteratorShifting & 1);
            *iteratorShifting = (*iteratorShifting + carryCurrent * basisCalc) >> 1;
        }
        charBin = carryNext ? '1' : '0';
        strBin.insert(strBin.begin(), 1, charBin);
    }
    return strBin;
}



} // namespace bigint



