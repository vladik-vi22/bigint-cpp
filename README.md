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

**Arithmetic:** `+`, `-`, `*`, `/`, `%`, `pow`, `log2`

**Modular arithmetic:** `powmod`, `inversemod`, `congruencemod`, `isCoprime`, `symbolJacobi`

**Bitwise:** `~`, `&`, `|`, `^`, `<<`, `>>`, circular shifts

**Number theory:** `gcd` (binary algorithm), `lcm`, `abs`

**I/O:** Construct from decimal/hex/binary strings, `std::vector<uint8_t>`, integers. Convert back to any format.

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
```

## Build

```bash
cmake -B build -DBIGINT_BUILD_TESTS=ON
cmake --build build --config Release
ctest --test-dir build -C Release
```

## Internals

- `std::vector<uint32_t>` storage, little-endian, base 2³²
- Schoolbook multiplication O(n²), binary GCD, square-and-multiply for powmod
- ~40 unit tests, Google Benchmark suite included

## License

MIT
