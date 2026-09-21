// libcccrt core - Shannon entropy
// Dependency-free (no Eigen, no heap). Suitable for embedded targets.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace cccrt {

// Shannon entropy (in bits) of a symbol sequence.
// Sorts `seq` in place and counts runs, so no histogram storage is needed.
template <typename Real = double, typename Sym>
Real shannonEntropyInPlace(Sym* seq, size_t n) {
    if (n == 0) return Real(0);
    std::sort(seq, seq + n);
    const Real scale = Real(1) / Real(n);
    Real H = 0;
    size_t i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && seq[j] == seq[i]) ++j;
        const Real p = Real(j - i) * scale;
        H -= p * std::log2(p);
        i = j;
    }
    return H;
}

// Non-destructive variant: `scratch` must hold at least n elements.
template <typename Real = double, typename Sym>
Real shannonEntropy(const Sym* seq, size_t n, Sym* scratch) {
    std::copy(seq, seq + n, scratch);
    return shannonEntropyInPlace<Real>(scratch, n);
}

} // namespace cccrt
