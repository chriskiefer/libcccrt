// PluginCccRPC.cpp
// Chris Kiefer (c.kiefer@sussex.ac.uk)

#include "SC_PlugIn.hpp"
#include "CccRPC.hpp"

static InterfaceTable* ft;

namespace Ccc {

CccRPC::CccRPC() {
    double maxWindowSize = 50.0; //make this a parameter
    mMaxWindowSize = static_cast<size_t>(maxWindowSize / 1000.0 * sampleRate());
    ringBuf.setSize(mMaxWindowSize);
    proj = RPC::createProjectionMatrix(static_cast<size_t>(in0(IN_HIGHDIM)), static_cast<size_t>(in0(IN_LOWDIM)));
    mCalcFunc = make_calc_function<CccRPC, &CccRPC::next>();
    next(1);
}

void CccRPC::next(int nSamples) {
    const float* input = in(IN_SIG);
    float* outbuf = out(0);
    double hopSize = 10.0; //make this a parameter
    double windowSize = 20.0; //make this a parameter

    for (int i = 0; i < nSamples; ++i) {
        using namespace std;
        ringBuf.push(input[i]);
        hopCounter++;
        size_t hopSizeInSamples = static_cast<size_t>(hopSize / 1000.0 * sampleRate());

        if (hopCounter >= min(hopSizeInSamples, mMaxWindowSize))
        {
            hopCounter = 0;
            size_t windowSizeInSamples =
                static_cast<size_t>(windowSize / 1000.0 * sampleRate());
            auto window = ringBuf.getBuffer(min(windowSizeInSamples, mMaxWindowSize));

            rpc = RPC::calc(proj, window, 10, 0.5);
        }
        outbuf[i] = rpc;
    }

}

} // namespace Ccc

PluginLoad(CccRPCUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<Ccc::CccRPC>(ft, "CccRPC", false);
}
