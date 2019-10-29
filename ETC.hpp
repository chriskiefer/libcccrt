#pragma once

#include <armadillo>
#include <iostream>
#include "shannonEntropy.hpp"
#include <tsl/hopscotch_map.h>

using namespace std;
using namespace arma;

static_assert(sizeof(sword) ==8); //findHFPair works on the assumption that arma::sword is 64 bit

struct ETC {

    union pair {
        struct {
            sword i1;
            sword i2;
        } __attribute__((packed));
        __int128 i128;

        bool operator==(const pair& other)
        {
            return i128 == other.i128;
        }
    };
    
    static ETC::pair makeETCPair (sword a, sword b) {ETC::pair p; p.i1 = a; p.i2 = b; return p;};

//    typedef unordered_map<__int128, unsigned int> pairFreqTable;
    typedef tsl::hopscotch_pg_map<__int128, unsigned int> pairFreqTable;

    static ivec symbolise(const vec &seq, const unsigned int bins)
    {
        return {0};
    }

    static ETC::pair findHFPair(const ivec &seq) {
        /*
         this implementation works very slightly differently from the matlab original where more that one pair wins the highest frequency
         - the matlab version depends on arbitrary behaviour of max(x)
         [m,indx]=max(Count_Array(:));
         i.e. chooses winner based on first position in 2d frequency matrix
         whereas this version chooses the first winner in the order of the array, which is slightly faster
         tests show this occasionally makes very minor differences in the results
         */
        ETC::pairFreqTable histo;
        ETC::pair winner;
        unsigned int highScore=0;

        unsigned int seqPos = 0;
        while (seqPos < seq.size()-1) {
            ETC::pair currPair = ETC::makeETCPair(seq[seqPos], seq[seqPos+1]);
            ETC::pairFreqTable::iterator it = histo.find(currPair.i128);
            unsigned int score;
            if (it == histo.end()) {
                histo.insert(std::make_pair(currPair.i128, 1));
                score=1;
            }else{
                score = it->second + 1;
                it.value() = score;
//                it->second = score;
            }
            if (score > highScore) {
                highScore = score;
                winner = currPair;
            }
            if (currPair.i1 == currPair.i2) {
                if (seqPos < seq.size()-1) {
                    if (seq[seqPos+2] == seq[seqPos]) {
                        seqPos++;
                    }
                }
            }
            seqPos++;
        }
        return winner;
    }
    

    static auto substitute(const ivec &seq, ETC::pair p) {
        sword replacementSymbol = max(seq) + 1;
        ivec newSeq = zeros<ivec>(seq.size());
        uword src=0, dest=0;
        unsigned int replaceCount=0;
        while(src < seq.size()) {
            if (src < seq.size()-1) {
                if (seq[src] == p.i1 && seq[src+1] == p.i2) {
                    newSeq[dest] = replacementSymbol;
                    src++;
                    replaceCount++;
                } else {
                    newSeq[dest] = seq[src];
                }
            }else{
                newSeq[dest] = seq[src];
            }
            src++;
            dest++;
        }
        // cout << dest << endl;
        newSeq = newSeq.rows(0,dest-1);
//        cout << newSeq << endl;
        return std::make_tuple(newSeq, replaceCount, replacementSymbol);
    }

    static double calcOld(const ivec &seq) {
        double N = 0; //ETC measure
        double Hnew = shannonEntropy::calc(seq);
        ivec newSeq = seq;
        
        //todo: this can be optimised by continually editing the shannon histo instead of redoing it every time
        while(Hnew >1e-6 && newSeq.size() > 1) {
            ETC::pair hfPair = ETC::findHFPair(newSeq);
            auto [newSeqRepl, replaceCount, replaceSym] = ETC::substitute(newSeq, hfPair);
            Hnew = shannonEntropy::calc(newSeqRepl);
            newSeq = newSeqRepl;
            N++;
            // cout << newSeq << endl;
            // cout << N << ", " << Hnew << endl;
        }
        
        return N;
    }
    
    static double calc(const ivec &seq) {
        double N = 0; //ETC measure
//        cout << seq << endl;
        if (seq.size() > 1) {
            shannonEntropy::histoMap histo = shannonEntropy::calcDistribution(seq);
            double Hnew = shannonEntropy::calcProbability(histo, seq);

            ivec newSeq = seq;
            
            while(Hnew >1e-6 && newSeq.size() > 1) {
//                cout << newSeq.size() << endl;
//                cout << newSeq.t() << endl;
                ETC::pair hfPair = ETC::findHFPair(newSeq);
                auto [newSeqRepl, replaceCount, replaceSym] = ETC::substitute(newSeq, hfPair);
    //            cout << newSeqRepl << endl;
                //reduce counts of replacement pair
                shannonEntropy::histoMap::iterator it = histo.find(hfPair.i1);
                it->second -= replaceCount;
                if (it->second == 0) {
                    //remove from the histo
                    histo.erase(it);
                }
                it = histo.find(hfPair.i2);
                it->second -= replaceCount;
                if (it->second == 0) {
                    //remove from the histo
                    histo.erase(it);
                }
                //add the new symbol into the histogram
                auto histoEntry = make_pair(replaceSym, replaceCount);
                histo.insert(histoEntry);
                
                
                Hnew = shannonEntropy::calcProbability(histo, newSeqRepl);
                newSeq = newSeqRepl;
                N++;
                // cout << newSeq << endl;
                // cout << N << ", " << Hnew << endl;
            }
            N /= (seq.size() -1);
        }
        return N;
    }
    
    static double calcJoint(const ivec& seq1, const ivec& seq2) {
        
        ivec combSeq;
        combSeq.set_size(seq1.size());
        for (uword i=0; i < seq1.size(); i++) {
            combSeq[i] =(sword)( seq1[i] | (seq2[i] << 32));
        }
        return ETC::calc(combSeq);
    }

};
