// libcccrt core - Sevcik fractal dimension
// Dependency-free (no Eigen, no heap). Suitable for embedded targets.
#pragma once

#include <cmath>
#include <cstddef>

namespace cccrt {

// Sevcik fractal dimension of a waveform.
// references:
//   https://arxiv.org/pdf/1003.5266.pdf
//   https://neuropsychology.github.io/NeuroKit/_modules/neurokit2/complexity/fractal_sevcik.html
// Returns 1 (the dimension of a line) for sequences shorter than 2 samples.
template <typename Real>
Real sevcik(const Real* x, size_t n) {
    if (n < 2) return Real(1);

    // scale 0 - 1
    Real mn = x[0], mx = x[0];
    for (size_t i = 1; i < n; ++i) {
        if (x[i] < mn) mn = x[i];
        if (x[i] > mx) mx = x[i];
    }
    const Real range = mx - mn;

    // path length of the normalised curve, with unit x-extent divided into n-1 steps
    Real dx = Real(1) / Real(n - 1);
    dx = dx * dx;
    Real L = 0;
    Real prev = (x[0] - mn) / range;
    for (size_t i = 1; i < n; ++i) {
        const Real cur = (x[i] - mn) / range;
        const Real dy = cur - prev;
        L += std::sqrt(dy * dy + dx);
        prev = cur;
    }
    return Real(1) + (std::log(L) / std::log(Real(2 * (n - 1))));
}

} // namespace cccrt
