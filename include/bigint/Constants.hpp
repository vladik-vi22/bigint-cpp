/**
 * @file Constants.hpp
 * @brief Public constants for the BigInt library.
 */

#pragma once

#include <cstdint>

namespace bigint {

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
inline constexpr uint8_t kDefaultOutputBase = kBaseDecimal;

/// @}

}  // namespace bigint
