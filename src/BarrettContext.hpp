/**
 * @file BarrettContext.hpp
 * @brief Barrett reduction context for fast modular reduction.
 * @internal This is an internal implementation header, not part of the public API.
 */

#ifndef BIGINT_BARRETT_CONTEXT_HPP
#define BIGINT_BARRETT_CONTEXT_HPP

#include <bigint/BigInt.hpp>

namespace bigint::internal {

/**
 * @brief Barrett reduction context for fast modular reduction.
 *
 * @details Barrett reduction replaces expensive division with multiplication
 *          by precomputing a scaled reciprocal of the modulus. This is useful
 *          when performing many modular operations with the same modulus.
 *
 *          Algorithm: To compute x mod n:
 *          1. Precompute mu = floor(2^(2k) / n) where k = bitLength(n)
 *          2. q = floor((x * mu) >> (2k))  (approximation of x / n)
 *          3. r = x - q * n
 *          4. While r >= n: r -= n  (at most 2 corrections needed)
 *
 *          Works for any modulus (unlike Montgomery which requires odd modulus).
 *
 * @see "Handbook of Applied Cryptography", Algorithm 14.42
 */
struct BarrettContext {
  /**
   * @brief Constructs Barrett context for given modulus.
   * @param modulus The modulus n (must be > 0).
   */
  explicit BarrettContext(const BigInt& modulus) : n_(modulus), k_(modulus.bitLength()) {
    // Compute mu = floor(2^(2k) / n)
    BigInt two_pow_2k = BigInt(1) << (2 * k_);
    mu_ = two_pow_2k / n_;
  }

  /**
   * @brief Computes x mod n using Barrett reduction.
   * @param x Value to reduce (must be non-negative and < n^2).
   * @return x mod n.
   */
  [[nodiscard]] BigInt reduce(const BigInt& x) const {
    // For small x, just use regular modulo
    if (x.bitLength() <= k_) {
      if (x < n_) {
        return x;
      }
      return x - n_;
    }

    // q = floor((x * mu) >> (2k))
    BigInt q = (x * mu_) >> (2 * k_);

    // r = x - q * n
    BigInt r = x - q * n_;

    // Correction: r may be slightly too large (at most 2n too large)
    while (r >= n_) {
      r -= n_;
    }

    // Handle negative case (shouldn't happen for valid input, but be safe)
    while (r < BigInt(0)) {
      r += n_;
    }

    return r;
  }

  /**
   * @brief Computes (a * b) mod n using Barrett reduction.
   * @param a First operand (should be < n).
   * @param b Second operand (should be < n).
   * @return (a * b) mod n.
   */
  [[nodiscard]] BigInt mulmod(const BigInt& a, const BigInt& b) const {
    return reduce(a * b);
  }

  BigInt n_;     ///< The modulus
  size_t k_;     ///< Bit length of modulus
  BigInt mu_;    ///< Precomputed reciprocal: floor(2^(2k) / n)
};

}  // namespace bigint::internal

#endif  // BIGINT_BARRETT_CONTEXT_HPP
