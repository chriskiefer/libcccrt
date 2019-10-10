#pragma once

#include <armadillo>
#include <iostream>
#include <unordered_map>
#include <cmath>
using namespace std;
using namespace arma;

struct shannonEntropy {
    
    typedef unordered_map<int, unsigned int> histoMap;
    typedef unordered_map<int, double> probMap;

    static shannonEntropy::histoMap calcDistribution(const ivec &seq) {
        shannonEntropy::histoMap histo;
        for (auto v: seq) {
            shannonEntropy::histoMap::iterator it = histo.find(v);
            if (it == histo.end()) {
                histo.insert(std::make_pair<int, unsigned int>(v, 1));
            }else{
                it->second = it->second + 1;
            }
        }
        return histo;
    }
    
    static double calcProbability(const shannonEntropy::histoMap &histo, const ivec &seq) {
        double scale = 1.0 / seq.size();
        double H=0;
        for(auto v: histo) {
            double prob = v.second * scale;
            H = H - (prob * log2(prob));
        }
        return H;
    }
    
    static double calc(const ivec &seq) {
        shannonEntropy::histoMap histo = shannonEntropy::calcDistribution(seq);
//        for(auto v: histo) {
//            cout << v.first << ": " << v.second << endl;
//        }
        double prob = shannonEntropy::calcProbability(histo, seq);
        return prob;
    }
};
