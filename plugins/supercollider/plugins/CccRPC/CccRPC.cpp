// PluginCccRPC.cpp
// Chris Kiefer (c.kiefer@sussex.ac.uk)

#include "SC_PlugIn.hpp"
#include "CccRPC.hpp"

static InterfaceTable* ft;

namespace Ccc {

CccRPC::CccRPC() {
    double maxWindowSize = in0(IN_MAXWINSIZE); 
    mMaxWindowSize = static_cast<size_t>(maxWindowSize / 1000.0 * sampleRate());
    ringBuf.setSize(mMaxWindowSize + 1);
    proj = RPC::createProjectionMatrix(static_cast<size_t>(in0(IN_HIGHDIM)), static_cast<size_t>(in0(IN_LOWDIM)));
    mCalcFunc = make_calc_function<CccRPC, &CccRPC::next>();
    next(1);
}

void CccRPC::next(int nSamples) {
    const float* input = in(IN_SIG);
    float* outbuf = out(0);
    double windowSize = in0(IN_WINSIZE); 
    double hopSize = in0(IN_HOPSIZE) * windowSize; 
    size_t windowSizeInSamples = static_cast<size_t>(windowSize / 1000.0 * sampleRate());
    size_t hopSizeInSamples = static_cast<size_t>(hopSize / 1000.0 * sampleRate());
    windowSize = sc_min(windowSizeInSamples, mMaxWindowSize);
    hopSizeInSamples = sc_min(hopSizeInSamples, mMaxWindowSize);

    for (int i = 0; i < nSamples; ++i) {
        ringBuf.push(input[i]);
        hopCounter++;

        if (hopCounter >= hopSizeInSamples)
        {
            hopCounter = 0;
            auto window = ringBuf.getBuffer(windowSize);

            rpc = RPC::calc(proj, window, in0(IN_RPCRES), in0(IN_RPCHOP));
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
