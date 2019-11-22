/*
 Kathpalia, A. and Nagaraj, N., 2019. Data-based intervention approach for Complexity-Causality measure. PeerJ Computer Science, 5, p.e196.
 */

#include "ETC.hpp"
#include <thread>
#include <future>


struct CCC {
    
    enum threading {MULTITHREAD, SINGLETHREAD};

    //DC(dx|xpast) = ETC(xpast + dx) - ETC(xpast)
    //equation 5 in the paper
    static std::tuple<double, bool> dynamicCC(const ivec &seq, size_t dx, size_t xpast, size_t step, threading threadType=CCC::MULTITHREAD) {
        size_t len = seq.size() -dx - xpast;
        auto calcCall = [&seq, dx, xpast, step, len]()  {
            double valCall=0;
            for(size_t i=0; i < len; i += step) {
                auto window = seq.subvec(i, i+ xpast + dx - 1);
                //                cout << "all: \n" << window.t() << endl;
                valCall += ETC::calc(window);
            }
            return valCall;
        };
        auto calcCpast = [&seq, xpast, step, len]() {
            double valCpast=0;
            for(size_t i=0; i < len; i += step) {
                auto window = seq.subvec(i, i+xpast-1);
                //                cout << "past: \n" << window.t() << endl;
                valCpast += ETC::calc(window);
            }
            return valCpast;
        };
        double Call=0;
        double Cpast=0;
        if (threadType == CCC::MULTITHREAD) {
            auto CallThread = std::async(std::launch::async, [&calcCall]() {
                return calcCall();
            });
            auto CpastThread = std::async(std::launch::async, [&calcCpast]() {
                return calcCpast();
            });
            Call = CallThread.get();
            Cpast = CpastThread.get();
        }else{
            Call = calcCall();
            Cpast = calcCpast();
        }
        int k = ceil(len/(double)step);
        double val = Call - Cpast;
        val = val / k;
        return std::make_tuple(val, Call > Cpast);
    }

    //equation 6 in the paper
    static double dynamicCCJoint(const ivec &X, const ivec &Y, size_t dx, size_t past, size_t step) {
        size_t len = X.size() -dx - past;
        int k=1;
        double val=0;
        for(size_t i=0; i < len; i += step) {
            auto CallThread = std::async(std::launch::async, [&X, &Y, dx, past, i]() {
                ivec window1 = X.subvec(i, i+ past + dx - 1);
                ivec window2 = X.subvec(i, i+ past + dx - 1);
                window2(span(0,past-1)) = Y(span(i,i+past-1));
//                auto window = seq1.subvec(i, i+ past + dx - 1);
//                cout << "all: \n" << window.t() << endl;
                return ETC::calcJoint(window1, window2);
            });
            auto CpastThread = std::async(std::launch::async, [&X, &Y, past, i]() {
                ivec window1 = X.subvec(i, i+past-1);
                ivec window2 = Y.subvec(i, i+past-1);
//                cout << "past: \n" << window.t() << endl;
                return ETC::calcJoint(window1, window2);
            });
            double Call = CallThread.get();
            double Cpast = CpastThread.get();
//            cout << Call << "," << Cpast << endl;
//            cout << "intCC: " << Call << ", " << Cpast << ", " << (Call - Cpast) <<endl;
            val += Call - Cpast;
            k++;
        }
        val = val / (k-1);
        return val;

    }

    //equation 8 in the paper
    static double CCCausality(const ivec &effectSeq, const ivec &causeSeq, size_t dx, size_t past, size_t step) {
        auto dynCCSeq1Thread = std::async(std::launch::async, [&effectSeq, past, dx, step]() {
            return CCC::dynamicCC(effectSeq, dx, past, step);
        });

        auto dynCCSeqJointThread = std::async(std::launch::async, [&effectSeq, &causeSeq, past, dx, step]() {
            return CCC::dynamicCCJoint(effectSeq, causeSeq, dx, past, step);
        });
        auto dynCCResult = dynCCSeq1Thread.get();
        double dynCC = get<0>(dynCCResult);
        double dynCCJoint = dynCCSeqJointThread.get();

//        cout << "CCC: " << dynCC << ", " << dynCCJoint << endl;

        return dynCC - dynCCJoint;
    }
};
