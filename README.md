# bigint-cpp

C++20 arbitrary precision integer library with number-theoretic functions.

## What is this?

A self-contained BigInt implementation extracted from a cryptography learning project. No dependencies beyond the standard library.

**When to use this:**
- Learning how arbitrary precision arithmetic works
- Small projects where you want zero dependencies
- Prototyping cryptographic algorithms

**When to use something else:**
- Production cryptography → use [GMP](https://gmplib.org/) (has everything, battle-tested)
- Need just big integers → [Boost.Multiprecision](https://www.boost.org/doc/libs/release/libs/multiprecision/) works great
- Note: Boost.Multiprecision lacks `inversemod` and Jacobi symbol; GMP has both

## Features

**Arithmetic:** `+`, `-`, `*`, `/`, `%`, `pow`, `sqrt`, `log2`

**Modular arithmetic:** `powmod`, `inversemod`, `congruencemod`, `isCoprime`, `symbolJacobi`

**Bitwise:** `~`, `&`, `|`, `^`, `<<`, `>>`, circular shifts, `testBit`

**Number theory:** `gcd` (Euclidean algorithm), `lcm`, `abs`, `isProbablePrime`, `nextPrime`, `randomPrime`

**Random:** `randomBits`, `randomBelow`

**I/O:** Construct from decimal/hex/binary strings, `std::vector<uint8_t>`, `std::span<T>`, integers. Convert back to any format.

## Installation

### Using CMake FetchContent (Recommended)

```cmake
include(FetchContent)

FetchContent_Declare(
    bigint-cpp
    GIT_REPOSITORY https://github.com/vladik-vi22/bigint-cpp.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(bigint-cpp)

target_link_libraries(your_target PRIVATE bigint::bigint)
```

### Manual Build

```bash
git clone https://github.com/vladik-vi22/bigint-cpp.git
cd bigint-cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

## Quick Example

```cpp
#include <bigint/BigInt.hpp>
using namespace bigint;

BigInt a("123456789012345678901234567890");
BigInt b("987654321098765432109876543210");

auto sum = a + b;
auto product = a * b;

// RSA-style modular exponentiation
BigInt base("12345"), exp("65537"), mod("1000000007");
auto encrypted = powmod(base, exp, mod);

// Modular inverse (for RSA decryption, etc.)
auto inv = inversemod(BigInt(17), BigInt(3120));  // 17^(-1) mod 3120

// Number theory
auto g = gcd(a, b);
auto jacobi = symbolJacobi(BigInt(1001), BigInt(9907));  // returns -1, 0, or 1

// Prime generation
auto prime = BigInt::randomPrime(256);  // 256-bit random prime
```

## Build

```bash
cmake -B build -DBIGINT_BUILD_TESTS=ON
cmake --build build --config Release
ctest --test-dir build -C Release
```

## Benchmarks

### vs Boost.Multiprecision

```bash
cmake -B build -DBIGINT_BUILD_BOOST_BENCHMARKS=ON
cmake --build build --config Release
./build/benchmarks/bigint_boost_benchmarks
```

| Operation | bigint-cpp | Boost | Ratio |
|-----------|------------|-------|-------|
| Add | 76 ns | 57 ns | ~1.3x slower |
| Multiply | 123 ns | 91 ns | ~1.4x slower |
| Multiply (Toom-3, 8K+ bits) | 89 μs | 26 μs | ~3.4x slower |
| Divide | 545 ns | 1.1 μs | **~2x faster!** |
| PowMod (e=65537) | 9.7 μs | 8.8 μs | ~1.1x slower |
| PowMod (large exp) | 558 μs | 1.2 ms | **~2x faster!** |
| GCD | 1.0 μs | 481 ns | ~2x slower |

### vs GMP (the gold standard)

```bash
# Requires GMP: vcpkg install gmp:x64-windows (or apt install libgmp-dev)
cmake -B build -DBIGINT_BUILD_GMP_BENCHMARKS=ON -DCMAKE_TOOLCHAIN_FILE=[vcpkg]/scripts/buildsystems/vcpkg.cmake
cmake --build build --config Release
./build/benchmarks/bigint_gmp_benchmarks
```

| Operation | bigint-cpp | GMP | Ratio |
|-----------|------------|-----|-------|
| Add | 79 ns | 10 ns | ~8x slower |
| Multiply | 119 ns | 16 ns | ~7x slower |
| Multiply (Toom-3, 8K+ bits) | 90 μs | 5.4 μs | ~17x slower |
| Divide | 545 ns | 130 ns | ~4x slower |
| PowMod (e=65537) | 10 μs | 437 ns | ~23x slower |
| PowMod (large exp) | 600 μs | 67 μs | ~9x slower |
| GCD | 1.1 μs | 155 ns | ~7x slower |

*Tested with ~200-600 bit numbers on MSVC 19.50, Release build.*

**Takeaway:** Division uses Knuth's Algorithm D - faster than Boost! For small exponents (RSA public key e=65537), we match Boost. For large exponents (RSA private key), our optimized Montgomery squaring beats Boost by 2x! GMP uses hand-tuned assembly - we're pure C++.

## Internals

- `std::vector<uint32_t>` storage, little-endian, base 2³²
- Schoolbook O(n²) multiplication for small numbers, Karatsuba O(n^1.585) for medium, Toom-Cook 3-way O(n^1.465) for very large (>8192 bits)
- **Knuth's Algorithm D** for division (TAOCP Vol 2, Section 4.3.1)
- Euclidean GCD (leverages fast division)
- PowMod: square-and-multiply with Montgomery CIOS + optimized squaring (large odd moduli) and Barrett reduction (even/medium moduli)
- Miller-Rabin primality testing with deterministic witnesses for small numbers
- 141 unit tests, Google Benchmark suite included

## License

MIT
