// cccrpc~ - Random Projection Complexity for Max/MSP
// Chris Kiefer (c.kiefer@sussex.ac.uk)
//
// A port of the CccRPC SuperCollider UGen (plugins/supercollider/plugins/CccRPC),
// built on the Eigen-free core (core/rpc.hpp).
//
// Parameter names follow the SuperCollider UGen; in the paper's notation
// h = highDim, l = lowDim, alpha = rpchop * highDim (samples), beta = res.
//
//   cccrpc~ [highDim] [lowDim] [maxWinSize] [maxLowDim]
//     highDim    : projection window length in samples, h (default 10) - fixed at creation
//     lowDim     : initial projection dimensions, l (default 2)
//     maxWinSize : maximum analysis window in ms (default 500)         - fixed at creation
//     maxLowDim  : upper limit for @lowdim (default 8, or lowDim if larger) - fixed at creation
//   attributes
//     @lowdim    : projection dimensions, l (1 .. maxLowDim)
//     @winsize   : analysis window in ms (default 25)
//     @hopsize   : analysis hop as a fraction of winsize (default 0.5)
//     @res       : histogram resolution per dimension (default 5)
//     @rpchop    : projection hop as a fraction of highDim (default 0.5)
//     @downsample: average this many input samples into one before analysis
//                  (default 1). N times cheaper; the window still spans the
//                  same time, but RPC then measures the smoothed signal.
//     @normalize : 0 = raw cell count (default)
//                  1 = divided by the maximum possible count, min(hops, res^l)
//                  2 = divided by the mean count for white noise with the
//                      current parameters (recalibrated when they change)
//
// Signal in. Left outlet: signal holding the most recent RPC value, updated
// once per analysis hop. Right outlet: the same value as a float, sent from the
// scheduler thread after each hop (at most once per signal vector).

// the Max SDK includes <windows.h>, whose min/max macros break std::min/std::max
#ifndef NOMINMAX
#define NOMINMAX
#endif

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
    const size_t maxLowDim;
    const double maxWinMs;

    double sampleRate = 0;
    size_t maxWindowSamples = 0;
    // maxLowDim x highDim, row-major. A projection into l dimensions uses the
    // first l rows, so @lowdim can change without touching the matrix.
    std::vector<double> matrix;
    std::vector<double> ringStorage;
    cccrt::RingBuffer<double> ring;
    std::vector<double> window;
    std::vector<double> projScratch;
    std::vector<uint64_t> cellScratch;
    size_t hopCounter = 0;
    double rpc = 0;      // raw cell count
    double output = 0;   // rpc after normalisation
    // downsampling accumulator
    double accum = 0;
    size_t accumCount = 0;
    // noise calibration (@normalize 2), computed on the main thread with its own scratch
    double noiseRef = 1;
    std::vector<double> calibWindow;
    std::vector<double> calibProj;
    std::vector<uint64_t> calibCells;

    CccRpcState(size_t hd, size_t maxLd, double mw) : highDim(hd), maxLowDim(maxLd), maxWinMs(mw) {
        matrix.resize(maxLowDim * highDim);
        cccrt::rpc::makeProjectionMatrix(matrix.data(), maxLowDim, highDim, 42);
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
        projScratch.assign(maxLowDim * maxHops, 0.0);
        cellScratch.assign(maxHops, 0);
        calibWindow.assign(maxWindowSamples, 0.0);
        calibProj.assign(maxLowDim * maxHops, 0.0);
        calibCells.assign(maxHops, 0);
        hopCounter = 0;
        rpc = 0;
        output = 0;
        accum = 0;
        accumCount = 0;
    }
};

typedef struct _cccrpc {
    t_pxobject ob;
    CccRpcState* state;
    void* out_float;
    void* clock;   // defers float output from the perform routine to the scheduler
    void* calibQelem; // runs noise calibration on the main thread after a parameter change
    // attributes
    long normalize;
    long lowdim;
    double winsize; // ms
    double hopsize; // fraction of winsize
    long res;
    double rpchop;
    long downsample;
} t_cccrpc;

static t_class* cccrpc_class = nullptr;

// the analysis parameters derived from the attributes, shared by the perform routine and calibration
struct AnalysisParams {
    size_t ds;
    size_t windowSamples;
    size_t hopSamples;
    size_t resolution;
    size_t lowDim;
    double rpcHop;
};

static AnalysisParams cccrpc_params(const t_cccrpc* x) {
    const CccRpcState& st = *x->state;
    AnalysisParams p;
    p.ds = static_cast<size_t>(std::max<long>(1, x->downsample));
    const double analysisRate = st.sampleRate / p.ds; // window/hop are in ms of input time
    const double winMs = x->winsize;
    const double hopMs = x->hopsize * winMs;
    p.windowSamples = static_cast<size_t>(winMs / 1000.0 * analysisRate);
    p.hopSamples = static_cast<size_t>(hopMs / 1000.0 * analysisRate);
    p.windowSamples = std::min(std::max<size_t>(1, p.windowSamples), st.maxWindowSamples);
    p.hopSamples = std::min(std::max<size_t>(1, p.hopSamples), st.maxWindowSamples);
    p.resolution = static_cast<size_t>(std::max<long>(1, x->res));
    p.lowDim = std::min(static_cast<size_t>(std::max<long>(1, x->lowdim)), st.maxLowDim);
    p.rpcHop = x->rpchop;
    return p;
}

void* cccrpc_new(t_symbol* s, long argc, t_atom* argv);
void cccrpc_calibrate(t_cccrpc* x);
t_max_err cccrpc_attr_set(t_cccrpc* x, t_object* attr, long argc, t_atom* argv);
void cccrpc_free(t_cccrpc* x);
void cccrpc_assist(t_cccrpc* x, void* b, long m, long a, char* s);
void cccrpc_tick(t_cccrpc* x);
void cccrpc_dsp64(t_cccrpc* x, t_object* dsp64, short* count, double samplerate, long maxvectorsize, long flags);
void cccrpc_perform64(t_cccrpc* x, t_object* dsp64, double** ins, long numins, double** outs, long numouts,
                      long sampleframes, long flags, void* userparam);

void C74_EXPORT ext_main(void* r) {
    t_class* c = class_new("cccrpc~", (method)cccrpc_new, (method)cccrpc_free, (long)sizeof(t_cccrpc), 0L, A_GIMME, 0);

    class_addmethod(c, (method)cccrpc_dsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method)cccrpc_assist, "assist", A_CANT, 0);

    CLASS_ATTR_LONG(c, "normalize", 0, t_cccrpc, normalize);
    CLASS_ATTR_FILTER_CLIP(c, "normalize", 0, 2);
    CLASS_ATTR_ENUMINDEX3(c, "normalize", 0, "Raw", "Maximum", "White Noise");
    CLASS_ATTR_LABEL(c, "normalize", 0, "Normalize Output");

    CLASS_ATTR_LONG(c, "lowdim", 0, t_cccrpc, lowdim);
    CLASS_ATTR_FILTER_MIN(c, "lowdim", 1);
    CLASS_ATTR_LABEL(c, "lowdim", 0, "Projection Dimensions (l)");

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

    CLASS_ATTR_LONG(c, "downsample", 0, t_cccrpc, downsample);
    CLASS_ATTR_FILTER_MIN(c, "downsample", 1);
    CLASS_ATTR_LABEL(c, "downsample", 0, "Downsample Factor");

    // every analysis parameter change re-runs the noise calibration
    for (const char* name : {"lowdim", "winsize", "hopsize", "res", "rpchop", "downsample"}) {
        CLASS_ATTR_ACCESSORS(c, name, NULL, cccrpc_attr_set);
    }

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
    long maxLowDim = 8;
    if (nPositional > 0) highDim = atom_getlong(argv);
    if (nPositional > 1) lowDim = atom_getlong(argv + 1);
    if (nPositional > 2) maxWinMs = atom_getfloat(argv + 2);
    if (nPositional > 3) maxLowDim = atom_getlong(argv + 3);
    highDim = std::max<long>(1, highDim);
    lowDim = std::max<long>(1, lowDim);
    if (maxWinMs <= 0.0) maxWinMs = 500.0;
    maxLowDim = std::max(std::max<long>(1, maxLowDim), lowDim);

    x->normalize = 0;
    x->lowdim = lowDim;
    x->winsize = 25.0;
    x->hopsize = 0.5;
    x->res = 5;
    x->rpchop = 0.5;
    x->downsample = 1;
    x->state = new CccRpcState(static_cast<size_t>(highDim), static_cast<size_t>(maxLowDim), maxWinMs);

    dsp_setup((t_pxobject*)x, 1);
    // outlets are created right to left
    x->out_float = floatout((t_object*)x);
    outlet_new((t_object*)x, "signal");
    x->clock = clock_new(x, (method)cccrpc_tick);
    x->calibQelem = qelem_new(x, (method)cccrpc_calibrate);
    attr_args_process(x, (short)argc, argv);
    return x;
}

void cccrpc_free(t_cccrpc* x) {
    dsp_free((t_pxobject*)x);
    object_free(x->clock);
    qelem_free(x->calibQelem);
    delete x->state;
    x->state = nullptr;
}

void cccrpc_assist(t_cccrpc* x, void* b, long m, long a, char* s) {
    if (m == ASSIST_INLET) {
        snprintf(s, 256, "(signal) Input");
    } else if (a == 0) {
        snprintf(s, 256, "(signal) Random projection complexity");
    } else {
        snprintf(s, 256, "(float) Random projection complexity, once per hop");
    }
}

void cccrpc_tick(t_cccrpc* x) {
    outlet_float(x->out_float, x->state->output);
}

// shared setter for the analysis attributes: store the value, then recalibrate
t_max_err cccrpc_attr_set(t_cccrpc* x, t_object* attr, long argc, t_atom* argv) {
    if (argc < 1 || !argv) return MAX_ERR_NONE;
    t_symbol* name = (t_symbol*)object_method(attr, gensym("getname"));
    if (name == gensym("lowdim")) x->lowdim = std::max<long>(1, atom_getlong(argv));
    else if (name == gensym("winsize")) x->winsize = std::max(0.0, atom_getfloat(argv));
    else if (name == gensym("hopsize")) x->hopsize = std::max(0.0, atom_getfloat(argv));
    else if (name == gensym("res")) x->res = std::max<long>(1, atom_getlong(argv));
    else if (name == gensym("rpchop")) x->rpchop = std::min(1.0, std::max(0.0, atom_getfloat(argv)));
    else if (name == gensym("downsample")) x->downsample = std::max<long>(1, atom_getlong(argv));
    if (x->calibQelem) qelem_set(x->calibQelem);
    return MAX_ERR_NONE;
}

// Reference value for @normalize 2: mean RPC of white noise windows with the
// current parameters. Runs on the main thread (qelem / dsp64) with its own scratch.
void cccrpc_calibrate(t_cccrpc* x) {
    CccRpcState& st = *x->state;
    if (st.sampleRate <= 0 || st.calibWindow.empty()) return; // dsp64 will calibrate once it knows the sample rate
    const AnalysisParams p = cccrpc_params(x);
    cccrt::rpc::Pcg32 rng(1234); // fixed seed: the same parameters always give the same reference
    const int nWindows = 8;
    double sum = 0;
    for (int w = 0; w < nWindows; ++w) {
        for (size_t i = 0; i < p.windowSamples; ++i) st.calibWindow[i] = rng.uniform01<double>() * 2.0 - 1.0;
        sum += cccrt::rpc::calc(st.matrix.data(), p.lowDim, st.highDim, st.calibWindow.data(), p.windowSamples,
                                p.resolution, p.rpcHop, st.calibProj.data(), st.calibCells.data());
    }
    st.noiseRef = std::max(1.0, sum / nWindows);
}

void cccrpc_dsp64(t_cccrpc* x, t_object* dsp64, short* count, double samplerate, long maxvectorsize, long flags) {
    x->state->prepare(samplerate);
    cccrpc_calibrate(x);
    object_method(dsp64, gensym("dsp_add64"), x, cccrpc_perform64, 0, NULL);
}

void cccrpc_perform64(t_cccrpc* x, t_object* dsp64, double** ins, long numins, double** outs, long numouts,
                      long sampleframes, long flags, void* userparam) {
    CccRpcState& st = *x->state;
    const double* in = ins[0];
    double* out = outs[0];

    // read the modulatable parameters once per block, as the UGen does
    const AnalysisParams p = cccrpc_params(x);
    const size_t ds = p.ds;
    const size_t windowSamples = p.windowSamples;
    const size_t hopSamples = p.hopSamples;
    const size_t resolution = p.resolution;
    const double rpcHop = p.rpcHop;
    const size_t lowDim = p.lowDim;
    const long normalize = x->normalize;
    bool newValue = false;

    for (long i = 0; i < sampleframes; ++i) {
        // average ds input samples into one analysis sample
        st.accum += in[i];
        if (++st.accumCount >= ds) {
            st.ring.push(st.accum / static_cast<double>(st.accumCount));
            st.accum = 0;
            st.accumCount = 0;
            if (++st.hopCounter >= hopSamples) {
                st.hopCounter = 0;
                st.ring.copyLatest(st.window.data(), windowSamples);
                st.rpc = cccrt::rpc::calc(st.matrix.data(), lowDim, st.highDim, st.window.data(), windowSamples,
                                          resolution, rpcHop, st.projScratch.data(), st.cellScratch.data());
                if (normalize == 1) {
                    const double mx = cccrt::rpc::maxOccupiedCells(windowSamples, st.highDim, rpcHop, resolution, lowDim);
                    st.output = mx > 0 ? st.rpc / mx : 0;
                } else if (normalize == 2) {
                    st.output = st.rpc / st.noiseRef;
                } else {
                    st.output = st.rpc;
                }
                newValue = true;
            }
        }
        out[i] = st.output;
    }
    if (newValue) clock_delay(x->clock, 0);
}
