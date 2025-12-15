/**
 * @file MontgomeryContext.cpp
 * @brief Implementation of Montgomery multiplication context.
 */

#include "MontgomeryContext.hpp"

namespace bigint::internal {

void MontgomeryContext::multiply(const WordVec& a, const WordVec& b, WordVec& result) const {
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

void MontgomeryContext::square(const WordVec& a, WordVec& result) const {
  std::vector<Word> t(2 * k_ + 1, 0);

  computeOffDiagonalTerms(a, t);
  doubleInPlace(t);
  addDiagonalTerms(a, t);
  montgomeryReduce(t);
  extractResult(t, result);
}

MontgomeryContext::Word MontgomeryContext::computeNegInverse(Word n0) noexcept {
  // Newton's method converges to n0^(-1) mod 2^32 in 5 iterations
  constexpr int kNewtonIterations = 5;
  Word inv = 1;
  for (int i = 0; i < kNewtonIterations; ++i) {
    inv *= 2 - n0 * inv;
  }
  return static_cast<Word>(-static_cast<int32_t>(inv));
}

BigInt MontgomeryContext::wordVecToBigInt(const WordVec& words) {
  // Reverse to big-endian for span constructor
  WordVec big_endian(words.rbegin(), words.rend());

  // Remove leading zeros
  auto first_nonzero = std::ranges::find_if(big_endian, [](Word w) { return w != 0; });
  if (first_nonzero != big_endian.begin()) {
    big_endian.erase(big_endian.begin(), first_nonzero);
  }
  if (big_endian.empty()) {
    big_endian.push_back(0);
  }

  return BigInt(std::span<const Word>(big_endian), true);
}

MontgomeryContext::WordVec MontgomeryContext::bigIntToWordVec(const BigInt& value, size_t target_size) {
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

MontgomeryContext::WordVec MontgomeryContext::computeRSquared() const {
  // R = 2^(32*k), so R^2 = 2^(64*k)
  WordVec r2(2 * k_ + 1, 0);
  r2[2 * k_] = 1;  // R^2 = 2^(64*k)
  return reduceModN(r2);
}

MontgomeryContext::WordVec MontgomeryContext::reduceModN(const WordVec& x) const {
  if (x.empty() || std::ranges::all_of(x, [](Word w) { return w == 0; })) {
    return WordVec(k_, 0);
  }

  BigInt big_x = wordVecToBigInt(x);
  BigInt big_n = wordVecToBigInt(n_);
  BigInt result = big_x % big_n;

  return bigIntToWordVec(result, k_);
}

void MontgomeryContext::conditionalSubtract(const WordVec& t, WordVec& result) const {
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

void MontgomeryContext::computeOffDiagonalTerms(const WordVec& a, WordVec& t) const {
  for (size_t i = 0; i < k_; ++i) {
    DWord carry = 0;
    for (size_t j = i + 1; j < k_; ++j) {
      DWord product = static_cast<DWord>(a[i]) * a[j] + t[i + j] + carry;
      t[i + j] = static_cast<Word>(product);
      carry = product >> kBitsPerDigit;
    }
    propagateCarry(t, i + k_, carry);
  }
}

void MontgomeryContext::doubleInPlace(WordVec& t) const {
  DWord carry = 0;
  for (size_t i = 0; i < t.size(); ++i) {
    DWord doubled = (static_cast<DWord>(t[i]) << 1) | carry;
    t[i] = static_cast<Word>(doubled);
    carry = doubled >> kBitsPerDigit;
  }
}

void MontgomeryContext::addDiagonalTerms(const WordVec& a, WordVec& t) const {
  DWord carry = 0;
  for (size_t i = 0; i < k_; ++i) {
    DWord product = static_cast<DWord>(a[i]) * a[i] + t[2 * i] + carry;
    t[2 * i] = static_cast<Word>(product);
    carry = product >> kBitsPerDigit;

    DWord sum = static_cast<DWord>(t[2 * i + 1]) + carry;
    t[2 * i + 1] = static_cast<Word>(sum);
    carry = sum >> kBitsPerDigit;
  }
  if (carry != 0) {
    t[2 * k_] += static_cast<Word>(carry);
  }
}

void MontgomeryContext::propagateCarry(WordVec& t, size_t start_idx, DWord carry) const {
  for (size_t j = start_idx; carry != 0 && j < t.size(); ++j) {
    DWord sum = static_cast<DWord>(t[j]) + carry;
    t[j] = static_cast<Word>(sum);
    carry = sum >> kBitsPerDigit;
  }
}

void MontgomeryContext::montgomeryReduce(WordVec& t) const {
  for (size_t i = 0; i < k_; ++i) {
    Word m = t[i] * n0_inv_;

    DWord carry = 0;
    for (size_t j = 0; j < k_; ++j) {
      DWord product = static_cast<DWord>(m) * n_[j] + t[i + j] + carry;
      t[i + j] = static_cast<Word>(product);
      carry = product >> kBitsPerDigit;
    }
    propagateCarry(t, i + k_, carry);
  }
}

void MontgomeryContext::extractResult(const WordVec& t, WordVec& result) const {
  WordVec t_reduced(k_ + 1, 0);
  for (size_t i = 0; i < k_; ++i) {
    t_reduced[i] = t[k_ + i];
  }
  t_reduced[k_] = t[2 * k_];
  conditionalSubtract(t_reduced, result);
}

}  // namespace bigint::internal
