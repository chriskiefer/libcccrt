// PluginCccRPC.hpp
// Chris Kiefer (c.kiefer@sussex.ac.uk)

#pragma once

#include "SC_PlugIn.hpp"
#include "../../../../RPC.hpp"
#include "../common/RingBuf.hpp"


namespace Ccc {

class CccRPC : public SCUnit {
public:
    CccRPC();

    // Destructor
    // ~CccRPC();

private:
    void next(int nSamples);

    enum Inputs { IN_SIG, IN_HIGHDIM, IN_LOWDIM, IN_RPCHOP, IN_RPCRES, IN_WINSIZE, IN_MAXWINSIZE, IN_HOPSIZE  };

    Eigen::MatrixXd proj;
    RingBuf<double> ringBuf;
    size_t           hopCounter = 0;
    double           rpc = 0;
    size_t           mMaxWindowSize;

};

} // namespace Ccc
