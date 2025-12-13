#pragma once

#include <cstdint>

namespace bigint {

// Forward declaration
class BigInt;

/// @name Base Constants
/// @{

/// Binary base for string conversion.
inline constexpr uint8_t kBaseBinary = 2;

/// Decimal base for string conversion.
inline constexpr uint8_t kBaseDecimal = 10;

/// Hexadecimal base for string conversion.
inline constexpr uint8_t kBaseHexadecimal = 16;

/// Default base for input string parsing.
inline constexpr uint8_t kDefaultInputBase = kBaseDecimal;

/// Default base for output string conversion.
inline constexpr int kDefaultOutputBase = kBaseDecimal;

/// @}

/// @name Numeric Constants
/// @brief Pre-defined BigInt constants for common values.
/// @{
namespace constants {

/// BigInt constant representing zero.
extern const BigInt kZero;

/// BigInt constant representing one.
extern const BigInt kOne;

/// BigInt constant representing two.
extern const BigInt kTwo;

/// BigInt constant representing three.
extern const BigInt kThree;

/// BigInt constant representing four.
extern const BigInt kFour;

/// BigInt constant representing five.
extern const BigInt kFive;

/// BigInt constant representing eight.
extern const BigInt kEight;

// Legacy names for backward compatibility (deprecated)
[[deprecated("Use kZero instead")]]
extern const BigInt& Zero;

[[deprecated("Use kOne instead")]]
extern const BigInt& One;

[[deprecated("Use kTwo instead")]]
extern const BigInt& Two;

[[deprecated("Use kThree instead")]]
extern const BigInt& Three;

[[deprecated("Use kFour instead")]]
extern const BigInt& Four;

[[deprecated("Use kFive instead")]]
extern const BigInt& Five;

[[deprecated("Use kEight instead")]]
extern const BigInt& Eight;

}  // namespace constants
/// @}

}  // namespace bigint

