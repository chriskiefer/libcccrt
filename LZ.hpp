#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "core/lz.hpp"

using ArrayXL = Eigen::Array<int64_t, Eigen::Dynamic, 1>; 

// Eigen-facing wrapper around core/lz.hpp
struct LZ {
    static size_t calc(const ArrayXL &seq) {
      std::vector<size_t> phraseStarts(seq.size() + 1);
      return cccrt::lempelZiv(seq.data(), seq.size(), phraseStarts.data());
    }

    /*Normalised 
    see Zhang, Y., Hao, J., Zhou, C., & Chang, K. (2009). Normalized Lempel-Ziv complexity and its application in bio-sequence analysis. Journal of Mathematical Chemistry, 46(4), 1203–1212. https://doi.org/10.1007/s10910-008-9512-2
    */
    static double calcNorm(const ArrayXL &seq) {
      std::vector<size_t> phraseStarts(seq.size() + 1);
      return cccrt::lempelZivNorm<double>(seq.data(), seq.size(), phraseStarts.data());
    }

};
