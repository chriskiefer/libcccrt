/*
 Kathpalia, A. and Nagaraj, N., 2019. Data-based intervention approach for Complexity-Causality measure. PeerJ Computer Science, 5, p.e196.
 */

#include "ETC.hpp"
#include <thread>
#include <future>
#include <Eigen/Dense>
 



struct CCC {
    
    enum threading {MULTITHREAD, SINGLETHREAD};

    //DC(dx|xpast) = ETC(xpast + dx) - ETC(xpast)
    //equation 5 in the paper
    static double dynamicCC(const ArrayXL &seq, size_t dx, size_t xpast, size_t step, threading threadType=CCC::MULTITHREAD) {
        size_t len = seq.size() -dx - xpast;
        auto calcCall = [&seq, dx, xpast, step, len]()  {
            double valCall=0;
            for(size_t i=0; i < len; i += step) {
                auto window = seq(Eigen::seq(i, i+ xpast + dx - 1));
                //                cout << "all: \n" << window.t() << endl;
                valCall += ETC::calc(window);
            }
            return valCall;
        };
        auto calcCpast = [&seq, xpast, step, len]() {
            double valCpast=0;
            for(size_t i=0; i < len; i += step) {
                ArrayXL window = seq(Eigen::seq(i, i+xpast-1));
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
        return val;
    }

    //equation 6 in the paper
    static double dynamicCCJoint(const ArrayXL &X, const ArrayXL &Y, size_t dx, size_t past, size_t step, threading threadType=CCC::MULTITHREAD) {
        size_t len = X.size() -dx - past;
        int k=1;
        double val=0;
        for(size_t i=0; i < len; i += step) {
            double Call, Cpast;
            //TODO: refactor repeated code using lamdas from multithread
            if (threadType==CCC::MULTITHREAD) {
                auto CallThread = std::async(std::launch::async, [&X, &Y, dx, past, i]() {
                    ArrayXL window1 = X(Eigen::seq(i, i+ past + dx - 1));
                    ArrayXL window2 = X(Eigen::seq(i, i+ past + dx - 1));
                    window2(Eigen::seq(0,past-1)) = Y(Eigen::seq(i,i+past-1));
                    return ETC::calcJoint(window1, window2);
                });
                auto CpastThread = std::async(std::launch::async, [&X, &Y, past, i]() {
                    ArrayXL window1 = X(Eigen::seq(i, i+past-1));
                    ArrayXL window2 = Y(Eigen::seq(i, i+past-1));
                    return ETC::calcJoint(window1, window2);
                });
                Call = CallThread.get();
                Cpast = CpastThread.get();
            }else{
                ArrayXL window1All= X(Eigen::seq(i, i+ past + dx - 1));
                ArrayXL window2All = X(Eigen::seq(i, i+ past + dx - 1));
                window2All(Eigen::seq(0,past-1)) = Y(Eigen::seq(i,i+past-1));
                Call = ETC::calcJoint(window1All, window2All);
                ArrayXL window1Past = X(Eigen::seq(i, i+past-1));
                ArrayXL window2Past = Y(Eigen::seq(i, i+past-1));
                Cpast = ETC::calcJoint(window1Past, window2Past);
            }
            val += Call - Cpast;
            k++;
        }
        val = val / (k-1);
        return val;

    }

    //equation 8 in the paper
    static std::tuple<double, unsigned int> CCCausality(const ArrayXL &effectSeq, const ArrayXL &causeSeq, size_t dx, size_t past, size_t step) {
        auto dynCCSeq1Thread = std::async(std::launch::async, [&effectSeq, past, dx, step]() {
            return CCC::dynamicCC(effectSeq, dx, past, step);
        });

        auto dynCCSeqJointThread = std::async(std::launch::async, [&effectSeq, &causeSeq, past, dx, step]() {
            return CCC::dynamicCCJoint(effectSeq, causeSeq, dx, past, step);
        });
        auto dynCCResult = dynCCSeq1Thread.get();
        double dynCC = dynCCResult;
        double dynCCJoint = dynCCSeqJointThread.get();
        double CCC = dynCC - dynCCJoint;
//        cout << "CCC: " << dynCC << ", " << dynCCJoint << endl;
        //CCC mode - see table S1, related to polarity of the terms, and relative magnitude
        unsigned int CCCMode;
        if (dynCC < 0) {
            if (dynCCJoint < 0) {
                if (CCC < 0) {
                    CCCMode = 0;
                }else{
                    CCCMode = 1;
                }
                
            }else{
                CCCMode = 2;
            }
        }else{
            if (dynCCJoint < 0) {
                CCCMode = 3;
            }else{
                if (CCC < 0) {
                    CCCMode = 4;
                }else{
                    CCCMode = 5;
                }
            }
        }

        return std::make_tuple(CCC, CCCMode);;
    }
};
