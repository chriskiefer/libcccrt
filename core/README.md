# libcccrt core

Dependency-free implementations of the lightweight libcccrt measures, for use
on microcontrollers (e.g. RP2350) or anywhere Eigen is unwanted.

* `shannon.hpp` — Shannon entropy
* `lz.hpp` — Lempel-Ziv complexity (raw and normalised)
* `sevcik.hpp` — Sevcik fractal dimension
* `rpc.hpp` — Random Projection Complexity
* `cccrt.hpp` — includes all of the above

Effort To Compress and Compression-Complexity Causality are not included; they
remain Eigen-based in the repository root.

The headers in the repository root (`shannonEntropy.hpp`, `LZ.hpp`,
`fractal.hpp`, `RPC.hpp`) are thin Eigen wrappers over these, so the Python
module, SuperCollider UGen and gtest suite all exercise the same code.

## Design

* Requires only a C++17 compiler and `<cmath>`/`<algorithm>`. Builds with
  `-fno-exceptions -fno-rtti` and newlib-nano.
* Functions take `const T* data, size_t n`; nothing is allocated. Where an
  algorithm needs working memory the caller passes a scratch buffer.
* Everything is templated on the scalar type. Use `float` on Cortex-M33: its
  FPU is single precision and `double` is emulated in software.
* `rpc::Fixed<Real, MaxDims, MaxWindow, MaxHops>` bundles the projection
  matrix and scratch space in statically sized members, for use from an audio
  callback without any heap access.

## Usage

```cpp
#include "cccrt.hpp"

// once, at startup
static cccrt::rpc::Fixed<float, 2, 32, 64> rpc;   // 2 dims, 32-sample projection window, up to 64 hops
rpc.init(2, 32);

// per analysis window (e.g. 500 samples from a ring buffer)
float rpcValue = rpc.calc(window, 500, /*resolution*/ 5, /*hop*/ 0.5f);
float fd       = cccrt::sevcik(window, 500);

// symbolic measures work on any integer type
uint8_t syms[500];  /* quantise window into syms */
static uint32_t lzScratch[501];              // n + 1 entries
static uint8_t  symScratch[500];             // n entries
float lz = cccrt::lempelZivNorm<float>(syms, 500, lzScratch);
float h  = cccrt::shannonEntropy<float>(syms, 500, symScratch);
// or, if the symbols can be destroyed:
float h2 = cccrt::shannonEntropyInPlace<float>(syms, 500);
```

The scratch sizes for the lower-level `rpc::calc` are `nDim * nHops` reals and
`nHops` `uint64_t`s, where `nHops = rpc::numHops(n, windowSize,
rpc::hopSizeInSamples(windowSize, hop))`.

### Randomness

The projection matrix is drawn with a small PCG32 + Box-Muller generator, so
values differ from those produced by the previous EigenRand-based code (and
between `float` and `double` builds). RPC is a relative measure so this is
normally irrelevant; if you need bit-identical results across platforms,
generate the matrix once with `rpc::makeProjectionMatrix` and embed it as a
`const` array — `rpc::calc` accepts any row-major `nDim x windowSize` buffer.

### Cost

RPC, Sevcik and Shannon entropy are O(n) per window (RPC is
`nDim * windowSize` multiply-adds per hop). Lempel-Ziv is roughly O(n²) and
will be the expensive one at audio rates on a microcontroller.

## Tests

`tests/core_test.cpp` is a standalone test program (no gtest) using the same
expected values as the main suite:

```sh
g++ -std=c++17 -O2 -fno-exceptions -fno-rtti core/tests/core_test.cpp -o core_test && ./core_test
```

It is also built as the `cccrt_core_test` target by the top-level CMake.
`tests/check-cortex-m33.sh` cross-compiles it with `arm-none-eabi-g++` using
RP2350 flags to check the headers build for the target.
