#include <iostream>
#include "shannonEntropy.hpp"
#include "ETC.hpp"
#include "CCC.hpp"
#include "LZ.hpp"
#include "RPC.hpp"
#include "fractal.hpp"
#include <Eigen/Dense>
#include "maximilian.h"

#include "npy.hpp" //https://github.com/llohse/libnpy
 
using Eigen::ArrayXi;

using namespace std;

maxiSample samp;


int main(int argc, char **argv) {
    cout << "CCC library benchmarks\n";
    samp.load("../jungle.wav");
    Eigen::Map<Eigen::VectorXd> data(samp.amplitudes.data(), samp.amplitudes.size()); 

    //perf testing
    auto proj = RPC::createProjectionMatrix(32,1);
    // auto data = Eigen::VectorXd::Random(1000,1);
    // cout << win << endl;
    auto dataCopy = Eigen::VectorXd(data);
    dataCopy = dataCopy.array() - dataCopy.minCoeff();
    dataCopy = dataCopy / dataCopy.maxCoeff() * 256.0;
    // cout << dataCopy << endl;

    auto dataSym = ArrayXL(data.size());
    for(size_t i=0; i < data.size(); i++) {
        dataSym[i] = static_cast<long>(dataCopy[i]);
    }
    // cout << dataSym << endl;

    const size_t runs = 1000000;
    size_t winLen = 500;
    const size_t numframes = static_cast<size_t>(floor(data.size() / winLen))-1;
    cout << numframes << endl;
    auto winSym = dataSym(Eigen::seqN(0,winLen));
    clock_t t = clock();
    vector<double> singleRunTimes(runs);
    // auto win = data.segment(0*winLen,winLen);
    for(int i=0; i < runs; i++) {
        if (i % 1000 == 0) {
            cout << "Runs: " << i << endl;
        }
        const size_t framenum = i % numframes;
        Eigen::Ref<Eigen::VectorXd> win = data.segment(framenum*winLen,winLen);
        clock_t tonce = clock();
        
        RPC::calc(proj, win, 100, 1);
        // shannonEntropy::calc(winSym);
        // LZ::calc(winSym);
        // fractal::sevcik::calc(win);
        // ETC::calc(winSym);
        
        const double work_time_once = (clock() - tonce) / double(CLOCKS_PER_SEC) * 1000;
        singleRunTimes[i] = work_time_once;
    }
    const double work_time = (clock() - t) / double(CLOCKS_PER_SEC) * 1000;
    cout << work_time << " ms" << endl;

    const std::vector<long unsigned> shape{1, runs};
    npy::SaveArrayAsNumpy("RPC0.npy", false, shape.size(), shape.data(), singleRunTimes);
    return 0;
}

