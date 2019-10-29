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
        cout << "----------\n" << seq.t() << endl;
        for(size_t i=0; i < len; i += step) {
            auto CallThread = std::async(std::launch::async, [&seq, dx, xpast, i]() {
                auto window = seq.subvec(i, i+ xpast + dx - 1);
                cout << "all: \n" << window.t() << endl;
                return ETC::calc(window);
            });
            double Call = CallThread.get();
            auto CpastThread = std::async(std::launch::async, [&seq, xpast, i]() {
                auto window = seq.subvec(i, i+xpast-1);
                cout << "past: \n" << window.t() << endl;
                return ETC::calc(window);
            });
            double Cpast = CpastThread.get();
            cout << Call << "," << Cpast << endl;
            val += Call - Cpast;
            k++;
        }
        if (k > 1) val = val / (k-1);
        return val;
    }
};
