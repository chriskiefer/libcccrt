#pragma once

#include <iostream>
#include <Eigen/Dense>

using ArrayXL = Eigen::Array<int64_t, Eigen::Dynamic, 1>; 

//calc lempel ziv
struct LZ {
    static size_t calc(const ArrayXL &seq) {
      vector<ArrayXL> dictionary;
      size_t wordStart=0;
      size_t wordEnd=1;
      while (wordEnd <= seq.size()) {
        auto word = seq.segment(wordStart, wordEnd - wordStart);
        // //is it in the dictionary
        auto findResult = std::find_if(dictionary.begin(), dictionary.end(), [word] (const ArrayXL& s) { 
            bool eq = false;
            if (s.size() == word.size()) {
              eq=true;
              for(size_t i=0; i < s.size(); i++) {
                if (word[i] != s[i]) {
                  eq = false;
                  break;
                }
              }
            }
            return eq;
        });
        if (findResult == dictionary.end()) {
          dictionary.push_back(word);
          wordStart = wordEnd;
          wordEnd = wordStart+1;
        }
        else{
          wordEnd++;
        }
      }
      return dictionary.size();
    }

    /*Normalised 
    see Zhang, Y., Hao, J., Zhou, C., & Chang, K. (2009). Normalized Lempel-Ziv complexity and its application in bio-sequence analysis. Journal of Mathematical Chemistry, 46(4), 1203–1212. https://doi.org/10.1007/s10910-008-9512-2
    */
    static double calcNorm(const ArrayXL &seq) {
      const auto n = seq.size();
      return LZ::calc(seq) / (n / log(n));
    }

};
