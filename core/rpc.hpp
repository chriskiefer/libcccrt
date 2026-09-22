// libcccrt core - Random Projection Complexity
// Dependency-free (no Eigen, no heap). Suitable for embedded targets.
//
// Kiefer, C. 2023. Dynamical complexity measurement with random projection:
// a metric optimised for realtime signal processing. Sound and Music Computing.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace cccrt {
namespace rpc {

// Minimal PCG32 generator (O'Neill), used only to build the projection matrix.
// Values differ from the EigenRand mt19937 generator used previously; RPC is a
// relative measure so this only matters if you need bit-identical numbers
// across platforms - in that case generate the matrix once and embed it.
struct Pcg32 {
    uint64_t state;
    uint64_t inc;

    explicit Pcg32(uint64_t seed, uint64_t stream = 1)
        : state(0), inc((stream << 1u) | 1u) {
        next();
        state += seed;
        next();
    }

    uint32_t next() {
        const uint64_t old = state;
        state = old * 6364136223846793005ULL + inc;
        const uint32_t xorshifted = static_cast<uint32_t>(((old >> 18u) ^ old) >> 27u);
        const uint32_t rot = static_cast<uint32_t>(old >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
    }

    // uniform in [0, 1) with 24 bits of resolution (exact in float)
    template <typename Real>
    Real uniform01() {
        return Real(next() >> 8) * Real(1.0 / 16777216.0);
    }

    // standard normal via Box-Muller
    template <typename Real>
    Real normal() {
        Real u1;
        do {
            u1 = uniform01<Real>();
        } while (u1 <= Real(0));
        const Real u2 = uniform01<Real>();
        return std::sqrt(Real(-2) * std::log(u1)) * std::cos(Real(6.283185307179586) * u2);
    }
};

// Fill `out` (row-major, nDim x windowSize) with N(0, 1/nDim) samples.
template <typename Real>
void makeProjectionMatrix(Real* out, size_t nDim, size_t windowSize, uint64_t seed = 42) {
    Pcg32 rng(seed);
    const Real sd = std::sqrt(Real(1) / Real(nDim));
    for (size_t i = 0; i < nDim * windowSize; ++i) {
        out[i] = rng.normal<Real>() * sd;
    }
}

inline size_t hopSizeInSamples(size_t windowSize, double hop) {
    return std::max(static_cast<size_t>(1), static_cast<size_t>(windowSize * hop));
}

inline size_t numHops(size_t n, size_t windowSize, size_t hopSize) {
    return n < windowSize ? 0 : ((n - windowSize) / hopSize) + 1;
}

// Upper bound on the value calc() can return for these parameters: each hop
// occupies at most one cell, and there are resolution^nDim cells. Useful for
// normalising the output to [0, 1].
inline double maxOccupiedCells(size_t n, size_t windowSize, double hop, size_t resolution, size_t nDim) {
    const double nHops = static_cast<double>(numHops(n, windowSize, hopSizeInSamples(windowSize, hop)));
    const double nCells = std::pow(static_cast<double>(resolution), static_cast<double>(nDim));
    return std::min(nHops, nCells);
}

// Flatten a multidimensional histogram bin index.
template <typename Idx>
uint64_t flatIndex(const Idx* tuple, size_t dims, size_t bound) {
    uint64_t index = static_cast<uint64_t>(tuple[0]);
    for (size_t i = 1; i < dims; ++i) {
        index *= bound;
        index += static_cast<uint64_t>(tuple[i]);
    }
    return index;
}

// Random projection complexity.
//   projectionMatrix : row-major nDim x windowSize (see makeProjectionMatrix)
//   data, n          : the signal
//   resolution       : histogram bins per dimension
//   hop              : hop between windows, as a fraction of windowSize
//   projScratch      : nDim * numHops(...) elements
//   cellScratch      : numHops(...) elements
// Returns the number of occupied histogram cells (0 if n < windowSize).
template <typename Real>
Real calc(const Real* projectionMatrix, size_t nDim, size_t windowSize,
          const Real* data, size_t n, size_t resolution, Real hop,
          Real* projScratch, uint64_t* cellScratch) {
    const size_t hopSize = hopSizeInSamples(windowSize, hop);
    const size_t nHops = numHops(n, windowSize, hopSize);
    if (nHops == 0) return Real(0);

    // project each window: projScratch[d * nHops + h]
    for (size_t h = 0; h < nHops; ++h) {
        const Real* window = data + h * hopSize;
        for (size_t d = 0; d < nDim; ++d) {
            const Real* row = projectionMatrix + d * windowSize;
            Real acc = 0;
            for (size_t k = 0; k < windowSize; ++k) {
                acc += row[k] * window[k];
            }
            projScratch[d * nHops + h] = acc;
        }
    }

    // translate each dimension to histogram bin indexes, in place
    for (size_t d = 0; d < nDim; ++d) {
        Real* row = projScratch + d * nHops;
        Real mn = row[0];
        for (size_t h = 1; h < nHops; ++h) mn = std::min(mn, row[h]);
        Real mx = 0;
        for (size_t h = 0; h < nHops; ++h) {
            row[h] -= mn;
            mx = std::max(mx, row[h]);
        }
        const Real scale = mx > Real(0) ? Real(1) / (mx * Real(1.000001)) : Real(1);
        for (size_t h = 0; h < nHops; ++h) {
            row[h] = std::floor(row[h] * scale * Real(resolution));
        }
    }

    // count occupied cells of the multidimensional histogram
    for (size_t h = 0; h < nHops; ++h) {
        uint64_t index = static_cast<uint64_t>(projScratch[h]);
        for (size_t d = 1; d < nDim; ++d) {
            index *= resolution;
            index += static_cast<uint64_t>(projScratch[d * nHops + h]);
        }
        cellScratch[h] = index;
    }
    std::sort(cellScratch, cellScratch + nHops);
    const uint64_t* last = std::unique(cellScratch, cellScratch + nHops);
    return Real(last - cellScratch);
}

// Convenience wrapper with fixed-capacity, statically allocated storage.
// Suitable for audio callbacks on microcontrollers.
template <typename Real, size_t MaxDims, size_t MaxWindow, size_t MaxHops>
class Fixed {
public:
    // Returns false (and leaves the object unusable) if the sizes exceed capacity.
    bool init(size_t nDim, size_t windowSize, uint64_t seed = 42) {
        if (nDim > MaxDims || windowSize > MaxWindow || nDim == 0 || windowSize == 0) {
            nDim_ = windowSize_ = 0;
            return false;
        }
        nDim_ = nDim;
        windowSize_ = windowSize;
        makeProjectionMatrix(matrix_, nDim_, windowSize_, seed);
        return true;
    }

    // Returns 0 if not initialised or if the data would need more than MaxHops windows.
    Real calc(const Real* data, size_t n, size_t resolution, Real hop = Real(0.5)) {
        if (nDim_ == 0) return Real(0);
        const size_t hopSize = hopSizeInSamples(windowSize_, hop);
        if (numHops(n, windowSize_, hopSize) > MaxHops) return Real(0);
        return rpc::calc(matrix_, nDim_, windowSize_, data, n, resolution, hop, proj_, cells_);
    }

    size_t nDim() const { return nDim_; }
    size_t windowSize() const { return windowSize_; }
    const Real* matrix() const { return matrix_; }

private:
    Real matrix_[MaxDims * MaxWindow];
    Real proj_[MaxDims * MaxHops];
    uint64_t cells_[MaxHops];
    size_t nDim_ = 0;
    size_t windowSize_ = 0;
};

} // namespace rpc
} // namespace cccrt
