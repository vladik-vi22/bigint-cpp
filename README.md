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

**Bitwise:** `~`, `&`, `|`, `^`, `<<`, `>>`, circular shifts

**Number theory:** `gcd` (binary algorithm), `lcm`, `abs`, `isProbablePrime`, `nextPrime`, `randomPrime`

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
cmake -B build -DBIGINT_BUILD_COMPARISON_BENCHMARKS=ON
cmake --build build --config Release
./build/benchmarks/bigint_comparison_benchmarks
```

| Operation | bigint-cpp | Boost | Ratio |
|-----------|------------|-------|-------|
| Add | 76 ns | 59 ns | ~1.3x slower |
| Multiply | 117 ns | 90 ns | ~1.3x slower |
| Divide | 67 μs | 1 μs | ~67x slower |
| PowMod | 462 μs | 9 μs | ~51x slower |
| GCD | 9 μs | 0.5 μs | ~18x slower |

### vs GMP (the gold standard)

```bash
# Requires GMP: vcpkg install gmp:x64-windows (or apt install libgmp-dev)
cmake -B build -DBIGINT_BUILD_GMP_BENCHMARKS=ON -DCMAKE_TOOLCHAIN_FILE=[vcpkg]/scripts/buildsystems/vcpkg.cmake
cmake --build build --config Release
./build/benchmarks/bigint_gmp_benchmarks
```

| Operation | bigint-cpp | GMP | Ratio |
|-----------|------------|-----|-------|
| Add | 77 ns | 10 ns | ~8x slower |
| Multiply | 120 ns | 17 ns | ~7x slower |
| Divide | 63 μs | 130 ns | ~480x slower |
| PowMod | 453 μs | 487 ns | ~930x slower |
| GCD | 8.8 μs | 158 ns | ~56x slower |

*Tested with ~600-bit numbers on MSVC 19.50, Release build.*

**Takeaway:** Addition/multiplication are competitive. Division and modular exponentiation have room for optimization (Newton-Raphson division, Montgomery multiplication). GMP uses hand-tuned assembly - we're pure C++.

## Internals

- `std::vector<uint32_t>` storage, little-endian, base 2³²
- Schoolbook O(n²) for small numbers, Karatsuba for large (threshold: 32 words)
- Binary GCD (Stein's algorithm), square-and-multiply for powmod
- Miller-Rabin primality testing with deterministic witnesses for small numbers
- 122 unit tests, Google Benchmark suite included

## License

MIT
