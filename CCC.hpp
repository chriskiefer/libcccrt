#include "ETC.hpp"
#include <thread>
#include <future>

//https://www.bogotobogo.com/cplusplus/C11/3_C11_Threading_Lambda_Functions.php

struct CCC {
    
    //DC(dx|xpast) = ETC(xpast + dx) - ETC(xpast)
    static double dynamicCC(const ivec &seq, size_t dx, size_t xpast, size_t step) {
        size_t len = seq.size() -dx - xpast;
        int k=0;
        double val=0;
//        cout << "----------\n" << seq.t() << endl;
        for(size_t i=0; i < len; i += step) {
            auto CallThread = std::async(std::launch::async, [&seq, dx, xpast, i]() {
                auto window = seq.subvec(i, i+ xpast + dx - 1);
//                cout << "all: \n" << window.t() << endl;
                return ETC::calc(window);
            });
            auto CpastThread = std::async(std::launch::async, [&seq, xpast, i]() {
                auto window = seq.subvec(i, i+xpast-1);
//                cout << "past: \n" << window.t() << endl;
                return ETC::calc(window);
            });
            double Call = CallThread.get();
            double Cpast = CpastThread.get();
            cout << Call << "," << Cpast << endl;
            val += Call - Cpast;
            k++;
        }
        if (k > 1) val = val / (k-1);
        return val;
    }
    
    static double dynamicCCJoint(const ivec &seq1, const ivec &seq2, size_t dx, size_t past, size_t step) {
        size_t len = seq1.size() -dx - past;
        int k=0;
        double val=0;
        for(size_t i=0; i < len; i += step) {
            auto CallThread = std::async(std::launch::async, [&seq1, &seq2, dx, past, i]() {
                ivec window1 = seq1.subvec(i, i+ past + dx - 1);
                ivec window2 = seq1.subvec(i, i+ past + dx - 1);
                window2(span(0,past-1)) = seq2(span(i,i+past-1));
//                auto window = seq1.subvec(i, i+ past + dx - 1);
//                cout << "all: \n" << window.t() << endl;
                return ETC::calcJoint(window1, window2);
            });
            double Call = CallThread.get();
            auto CpastThread = std::async(std::launch::async, [&seq1, &seq2, past, i]() {
                ivec window1 = seq1.subvec(i, i+past-1);
                ivec window2 = seq2.subvec(i, i+past-1);
//                cout << "past: \n" << window.t() << endl;
                return ETC::calcJoint(window1, window2);
            });
            double Cpast = CpastThread.get();
//            cout << Call << "," << Cpast << endl;
            val += Call - Cpast;
            k++;
        }
        if (k > 1) val = val / (k-1);
        return val;

    }
};
