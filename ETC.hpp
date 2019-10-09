#include <armadillo>
#include <iostream>
#include "shannonEntropy.hpp"

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

    typedef map<__int128, unsigned int> pairFreqTable; //note: unordered_map will probably run faster - does the order of pairs need to be deterministic?

    static ivec symbolise(const vec &seq, const unsigned int bins)
    {
        return {0};
    }


    //todo: to optimise, multi thread this?
    static ETC::pair findHFPair(const ivec &seq) {

        ETC::pairFreqTable histo;
        unsigned int seqPos = 0;
        while (seqPos < seq.size()-1) {
            ETC::pair currPair = ETC::makeETCPair(seq[seqPos], seq[seqPos+1]);
            ETC::pairFreqTable::iterator it = histo.find(currPair.i128);
            if (it == histo.end()) {
                histo.insert(std::make_pair(currPair.i128, 1));
            }else{
                it->second = it->second + 1;
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
        //find the highest freq pair in the frequency table
        ETC::pair winner;
        unsigned int highScore=0;
//        cout << "res: " << endl;
        for(auto v: histo) {
            ETC::pair p;
            p.i128 = v.first;
//            cout << p.i1 << ", " << p.i2 << ": \t" << v.second << endl;
            if (v.second > highScore) {
                highScore = v.second;
                winner = p;
            }
        }
        return winner;
    }

    static ivec substitute(const ivec &seq, ETC::pair p) {
        int replacementSymbol = max(seq) + 1;
        ivec newSeq = zeros<ivec>(seq.size());
        uword src=0, dest=0;
        while(src < seq.size()) {
            if (src < seq.size()-1) {
                if (seq[src] == p.i1 && seq[src+1] == p.i2) {
                    newSeq[dest] = replacementSymbol;
                    src++;
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
        return newSeq;
    }

    static double calc(const ivec &seq) {
        double N = 0; //ETC measure
        double Hnew = shannonEntropy::calc(seq);
        ivec newSeq = seq;

        //todo: this can be optimised by continually editing the shannon histo instead of redoing it every time
        while(Hnew >1e-6 && newSeq.size() > 1) {
            ETC::pair hfPair = ETC::findHFPair(newSeq);
            newSeq = ETC::substitute(newSeq, hfPair);
            Hnew = shannonEntropy::calc(newSeq);
            N++;
            // cout << newSeq << endl;
            // cout << N << ", " << Hnew << endl;
        }

        return N;
    }
};
