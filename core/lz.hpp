// libcccrt core - Lempel-Ziv complexity
// Dependency-free (no Eigen, no heap). Suitable for embedded targets.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace cccrt {

// Lempel-Ziv complexity: the number of distinct phrases the sequence parses into.
//
// The parsed phrases are contiguous in `seq`, so the dictionary is stored as
// phrase start offsets rather than copies of the words.
// `phraseStarts` must hold at least n + 1 entries.
template <typename Sym, typename Index = uint32_t>
size_t lempelZiv(const Sym* seq, size_t n, Index* phraseStarts) {
    size_t nPhrases = 0;
    size_t wordStart = 0;
    size_t wordEnd = 1;
    phraseStarts[0] = 0;
    while (wordEnd <= n) {
        const size_t len = wordEnd - wordStart;
        bool found = false;
        for (size_t k = 0; k < nPhrases && !found; ++k) {
            const size_t pStart = phraseStarts[k];
            const size_t pLen = phraseStarts[k + 1] - pStart;
            if (pLen == len) {
                found = std::equal(seq + wordStart, seq + wordEnd, seq + pStart);
            }
        }
        if (!found) {
            ++nPhrases;
            phraseStarts[nPhrases] = static_cast<Index>(wordEnd);
            wordStart = wordEnd;
            wordEnd = wordStart + 1;
        } else {
            ++wordEnd;
        }
    }
    return nPhrases;
}

// Normalised Lempel-Ziv complexity: LZ / (n / log(n)).
// Zhang, Y., Hao, J., Zhou, C., & Chang, K. (2009). Normalized Lempel-Ziv
// complexity and its application in bio-sequence analysis.
// Journal of Mathematical Chemistry, 46(4), 1203-1212.
template <typename Real = double, typename Sym, typename Index = uint32_t>
Real lempelZivNorm(const Sym* seq, size_t n, Index* phraseStarts) {
    const Real nr = Real(n);
    return Real(lempelZiv(seq, n, phraseStarts)) / (nr / std::log(nr));
}

} // namespace cccrt
