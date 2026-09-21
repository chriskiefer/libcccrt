#pragma once

#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "core/shannon.hpp"
using namespace std;
using ArrayXL = Eigen::Array<int64_t, Eigen::Dynamic, 1>; 


// Eigen-facing wrapper around core/shannon.hpp.
// The histogram-based functions are kept for ETC, which updates the histogram incrementally.
struct shannonEntropy {
    
    typedef unordered_map<int64_t, size_t> histoMap;


    static shannonEntropy::histoMap calcDistribution(const ArrayXL &seq) {
        shannonEntropy::histoMap histo;
        for (const uint64_t v: seq) {
            shannonEntropy::histoMap::iterator it = histo.find(v);
            if (it == histo.end()) {
                histo.insert(std::make_pair(v, 1));
            }else{
                it->second = it->second + 1;
            }
        }
        return histo;
    }
    
    static double calcProbability(const shannonEntropy::histoMap &histo, const ArrayXL &seq) {
        double scale = 1.0 / seq.size();
        double H=0;
        for(auto v: histo) {
            double prob = v.second * scale;
            H = H - (prob * log2(prob));
        }
        return H;
    }
    
    static double calc(const ArrayXL &seq) {
        std::vector<int64_t> scratch(seq.size());
        return cccrt::shannonEntropy<double>(seq.data(), seq.size(), scratch.data());
    }
};
