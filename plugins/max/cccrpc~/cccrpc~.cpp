// cccrpc~ - Random Projection Complexity for Max/MSP
// Chris Kiefer (c.kiefer@sussex.ac.uk)
//
// A port of the CccRPC SuperCollider UGen (plugins/supercollider/plugins/CccRPC),
// built on the Eigen-free core (core/rpc.hpp).
//
//   cccrpc~ [highDim] [lowDim] [maxWinSize]
//     highDim    : projection window length in samples (default 10) - fixed at creation
//     lowDim     : projection dimensions (default 2)                - fixed at creation
//     maxWinSize : maximum analysis window in ms (default 500)      - fixed at creation
//   attributes
//     @winsize   : analysis window in ms (default 25)
//     @hopsize   : analysis hop as a fraction of winsize (default 0.5)
//     @res       : histogram resolution per dimension (default 5)
//     @rpchop    : projection hop as a fraction of highDim (default 0.5)
//
// Signal in, signal out. The output holds the most recent RPC value and is
// updated once per analysis hop.

#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"

#include "rpc.hpp"
#include "ringbuf.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>

// C++ state lives behind a pointer because Max allocates objects with object_alloc (no constructors).
struct CccRpcState {
    const size_t highDim;
    const size_t lowDim;
    const double maxWinMs;

    double sampleRate = 0;
    size_t maxWindowSamples = 0;
    std::vector<double> matrix; // lowDim x highDim, row-major
    std::vector<double> ringStorage;
    cccrt::RingBuffer<double> ring;
    std::vector<double> window;
    std::vector<double> projScratch;
    std::vector<uint64_t> cellScratch;
    size_t hopCounter = 0;
    double rpc = 0;

    CccRpcState(size_t hd, size_t ld, double mw) : highDim(hd), lowDim(ld), maxWinMs(mw) {
        matrix.resize(lowDim * highDim);
        cccrt::rpc::makeProjectionMatrix(matrix.data(), lowDim, highDim, 42);
    }

    // Called from dsp64 (main thread), so allocation is fine here.
    void prepare(double sr) {
        if (sr == sampleRate && maxWindowSamples > 0) return;
        sampleRate = sr;
        maxWindowSamples = std::max<size_t>(1, static_cast<size_t>(maxWinMs / 1000.0 * sr));
        ringStorage.assign(maxWindowSamples + 1, 0.0);
        ring.init(ringStorage.data(), ringStorage.size());
        window.assign(maxWindowSamples, 0.0);
        // scratch for the worst case: a projection hop of one sample over the largest window
        const size_t maxHops = std::max<size_t>(1, cccrt::rpc::numHops(maxWindowSamples, highDim, 1));
        projScratch.assign(lowDim * maxHops, 0.0);
        cellScratch.assign(maxHops, 0);
        hopCounter = 0;
        rpc = 0;
    }
};

typedef struct _cccrpc {
    t_pxobject ob;
    CccRpcState* state;
    // attributes
    double winsize; // ms
    double hopsize; // fraction of winsize
    long res;
    double rpchop;
} t_cccrpc;

static t_class* cccrpc_class = nullptr;

void* cccrpc_new(t_symbol* s, long argc, t_atom* argv);
void cccrpc_free(t_cccrpc* x);
void cccrpc_assist(t_cccrpc* x, void* b, long m, long a, char* s);
void cccrpc_dsp64(t_cccrpc* x, t_object* dsp64, short* count, double samplerate, long maxvectorsize, long flags);
void cccrpc_perform64(t_cccrpc* x, t_object* dsp64, double** ins, long numins, double** outs, long numouts,
                      long sampleframes, long flags, void* userparam);

void C74_EXPORT ext_main(void* r) {
    t_class* c = class_new("cccrpc~", (method)cccrpc_new, (method)cccrpc_free, (long)sizeof(t_cccrpc), 0L, A_GIMME, 0);

    class_addmethod(c, (method)cccrpc_dsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method)cccrpc_assist, "assist", A_CANT, 0);

    CLASS_ATTR_DOUBLE(c, "winsize", 0, t_cccrpc, winsize);
    CLASS_ATTR_FILTER_MIN(c, "winsize", 0.0);
    CLASS_ATTR_LABEL(c, "winsize", 0, "Window Size (ms)");

    CLASS_ATTR_DOUBLE(c, "hopsize", 0, t_cccrpc, hopsize);
    CLASS_ATTR_FILTER_MIN(c, "hopsize", 0.0);
    CLASS_ATTR_LABEL(c, "hopsize", 0, "Hop Size (fraction of window)");

    CLASS_ATTR_LONG(c, "res", 0, t_cccrpc, res);
    CLASS_ATTR_FILTER_MIN(c, "res", 1);
    CLASS_ATTR_LABEL(c, "res", 0, "Histogram Resolution");

    CLASS_ATTR_DOUBLE(c, "rpchop", 0, t_cccrpc, rpchop);
    CLASS_ATTR_FILTER_CLIP(c, "rpchop", 0.0, 1.0);
    CLASS_ATTR_LABEL(c, "rpchop", 0, "Projection Hop (fraction of highDim)");

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    cccrpc_class = c;
}

void* cccrpc_new(t_symbol* s, long argc, t_atom* argv) {
    t_cccrpc* x = (t_cccrpc*)object_alloc(cccrpc_class);
    if (!x) return nullptr;

    // positional args, then attributes
    const long nPositional = attr_args_offset((short)argc, argv);
    long highDim = 10;
    long lowDim = 2;
    double maxWinMs = 500.0;
    if (nPositional > 0) highDim = atom_getlong(argv);
    if (nPositional > 1) lowDim = atom_getlong(argv + 1);
    if (nPositional > 2) maxWinMs = atom_getfloat(argv + 2);
    highDim = std::max<long>(1, highDim);
    lowDim = std::max<long>(1, lowDim);
    if (maxWinMs <= 0.0) maxWinMs = 500.0;

    x->winsize = 25.0;
    x->hopsize = 0.5;
    x->res = 5;
    x->rpchop = 0.5;
    x->state = new CccRpcState(static_cast<size_t>(highDim), static_cast<size_t>(lowDim), maxWinMs);

    dsp_setup((t_pxobject*)x, 1);
    outlet_new((t_object*)x, "signal");
    attr_args_process(x, (short)argc, argv);
    return x;
}

void cccrpc_free(t_cccrpc* x) {
    dsp_free((t_pxobject*)x);
    delete x->state;
    x->state = nullptr;
}

void cccrpc_assist(t_cccrpc* x, void* b, long m, long a, char* s) {
    if (m == ASSIST_INLET) {
        snprintf(s, 256, "(signal) Input");
    } else {
        snprintf(s, 256, "(signal) Random projection complexity");
    }
}

void cccrpc_dsp64(t_cccrpc* x, t_object* dsp64, short* count, double samplerate, long maxvectorsize, long flags) {
    x->state->prepare(samplerate);
    object_method(dsp64, gensym("dsp_add64"), x, cccrpc_perform64, 0, NULL);
}

void cccrpc_perform64(t_cccrpc* x, t_object* dsp64, double** ins, long numins, double** outs, long numouts,
                      long sampleframes, long flags, void* userparam) {
    CccRpcState& st = *x->state;
    const double* in = ins[0];
    double* out = outs[0];

    // read the modulatable parameters once per block, as the UGen does
    const double winMs = x->winsize;
    const double hopMs = x->hopsize * winMs;
    size_t windowSamples = static_cast<size_t>(winMs / 1000.0 * st.sampleRate);
    size_t hopSamples = static_cast<size_t>(hopMs / 1000.0 * st.sampleRate);
    windowSamples = std::min(std::max<size_t>(1, windowSamples), st.maxWindowSamples);
    hopSamples = std::min(std::max<size_t>(1, hopSamples), st.maxWindowSamples);
    const size_t resolution = static_cast<size_t>(std::max<long>(1, x->res));
    const double rpcHop = x->rpchop;

    for (long i = 0; i < sampleframes; ++i) {
        st.ring.push(in[i]);
        if (++st.hopCounter >= hopSamples) {
            st.hopCounter = 0;
            st.ring.copyLatest(st.window.data(), windowSamples);
            st.rpc = cccrt::rpc::calc(st.matrix.data(), st.lowDim, st.highDim, st.window.data(), windowSamples,
                                      resolution, rpcHop, st.projScratch.data(), st.cellScratch.data());
        }
        out[i] = st.rpc;
    }
}
