/**
 * @file MontgomeryContext.hpp
 * @brief Montgomery multiplication context for modular exponentiation.
 * @internal This is an internal implementation header, not part of the public API.
 */

#ifndef BIGINT_MONTGOMERY_CONTEXT_HPP
#define BIGINT_MONTGOMERY_CONTEXT_HPP

#include <bigint/BigInt.hpp>

#include <algorithm>
#include <cstdint>
#include <ranges>
#include <vector>

namespace bigint {
namespace internal {

/// Number of bits per digit (32-bit words)
constexpr size_t kBitsPerDigit = 32;

/**
 * @brief Montgomery multiplication context for modular exponentiation.
 *
 * @details Implements the CIOS (Coarsely Integrated Operand Scanning) algorithm
 * which combines multiplication and reduction in a single pass, avoiding
 * temporary BigInt allocations.
 *
 * Montgomery form: x~ = x * R mod N, where R = 2^(32*k) and k = word count of N.
 * Multiplication in Montgomery form: montMul(a~, b~) = a*b*R mod N
 *
 * @note This struct is used exclusively by the powmod function for efficient
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
   * @brief Converts a value to Montgomery form: x~ = x * R mod N.
   * @param x Input value (must be < N, zero-padded to k words).
   * @return x in Montgomery form.
   */
  [[nodiscard]] WordVec toMontgomery(const WordVec& x) const {
    WordVec result(k_);
    multiply(x, r_squared_, result);
    return result;
  }

  /**
   * @brief Converts from Montgomery form back to normal: x = x~ * R^(-1) mod N.
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
  void multiply(const WordVec& a, const WordVec& b, WordVec& result) const;

  /**
   * @brief Montgomery squaring: result = a^2 * R^(-1) mod N.
   * @param a Operand in Montgomery form.
   * @param[out] result Square in Montgomery form.
   * @details Optimized squaring exploits symmetry: a[i]*a[j] = a[j]*a[i] for i!=j.
   */
  void square(const WordVec& a, WordVec& result) const;

  /**
   * @brief Computes -N^(-1) mod 2^32 using Newton's method.
   * @param n0 Least significant word of N (must be odd).
   * @return The value -N^(-1) mod 2^32.
   */
  [[nodiscard]] static Word computeNegInverse(Word n0) noexcept;

  /**
   * @brief Converts little-endian word vector to BigInt.
   * @param words Little-endian word vector.
   * @return BigInt representation (always positive).
   */
  [[nodiscard]] static BigInt wordVecToBigInt(const WordVec& words);

  /**
   * @brief Converts BigInt to little-endian word vector.
   * @param value BigInt to convert.
   * @param target_size Desired size of output vector (padded with zeros).
   * @return Little-endian word vector.
   */
  [[nodiscard]] static WordVec bigIntToWordVec(const BigInt& value, size_t target_size);

 private:
  [[nodiscard]] WordVec computeRSquared() const;
  [[nodiscard]] WordVec reduceModN(const WordVec& x) const;
  void conditionalSubtract(const WordVec& t, WordVec& result) const;
  void computeOffDiagonalTerms(const WordVec& a, WordVec& t) const;
  void doubleInPlace(WordVec& t) const;
  void addDiagonalTerms(const WordVec& a, WordVec& t) const;
  void propagateCarry(WordVec& t, size_t start_idx, DWord carry) const;
  void montgomeryReduce(WordVec& t) const;
  void extractResult(const WordVec& t, WordVec& result) const;

  WordVec n_;                    ///< Modulus N
  size_t k_;                     ///< Word count of N
  Word n0_inv_;                  ///< -N^(-1) mod 2^32
  WordVec r_squared_;            ///< R^2 mod N for conversion to Montgomery form
  mutable WordVec scratch_;      ///< Scratch buffer for CIOS algorithm
};

}  // namespace internal
}  // namespace bigint

#endif  // BIGINT_MONTGOMERY_CONTEXT_HPP
