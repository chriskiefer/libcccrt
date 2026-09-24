// cccrpc~ - Random Projection Complexity for Max/MSP
// Chris Kiefer (c.kiefer@sussex.ac.uk)
//
// A port of the CccRPC SuperCollider UGen (plugins/supercollider/plugins/CccRPC),
// built on the Eigen-free core (core/rpc.hpp).
//
// Parameter names follow the SuperCollider UGen; in the paper's notation
// h = highDim, l = lowDim, alpha = rpchop * highDim (samples), beta = res.
//
//   cccrpc~ [@attribute value ...]
//   attributes
//     @highdim   : projection window length in samples, h (default 10)
//     @lowdim    : projection dimensions, l (default 4)
//     @res       : histogram resolution per dimension (default 10)
//     @rpchop    : projection hop as a fraction of highdim (default 0.5)
//     @maxlowdim : @lowdim values up to this are preallocated (default 8); a
//                  larger @lowdim still works, but reallocates
//     @winsize   : analysis window in ms (default 25)
//     @hopsize   : analysis hop as a fraction of winsize (default 0.5)
//     @downsample: average this many input samples into one before analysis
//                  (default 1). N times cheaper; the window still spans the
//                  same time, but RPC then measures the smoothed signal.
//     @maxwinsize: longest @winsize allowed, in ms (default 500)
//     @normalize : 0 = raw cell count (default)
//                  1 = divided by the maximum possible count, min(hops, res^l)
//                  2 = divided by the mean count for white noise with the
//                      current parameters (recalibrated when they change)
//
//   Older patches may give cccrpc~ [highDim] [lowDim] [maxWinSize] [maxLowDim]
//   as arguments; these set the matching attributes, and a typed or saved
//   attribute takes precedence.
//
// Changing @highdim, @maxwinsize, or @lowdim beyond the allocated size rebuilds
// the analysis state on the main thread; the audio thread picks the new state
// up at its next signal vector, and the analysis restarts from an empty window.
//
// On res and lowdim: the output is capped at min(hops, res^lowdim) occupied
// cells. The defaults give 10^4 cells against a few hundred hops, so the hop
// count is the limit and the measure keeps discriminating at the top of its
// range. Lowering them (e.g. res 5, lowdim 2 = 25 cells) saturates on
// broadband material, where it shows up as the output going deaf to high
// frequencies - the value has already hit the ceiling.
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
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <new>
#include <vector>

// Attribute defaults. The STR() of each feeds CLASS_ATTR_DEFAULT_SAVE, so the
// inspector and cccrpc_new agree (Max does not apply attribute defaults to
// non-UI objects; cccrpc_new writes them into the struct).
#define DEFAULT_HIGHDIM 10
#define DEFAULT_LOWDIM 4
#define DEFAULT_RES 10
#define DEFAULT_RPCHOP 0.5
#define DEFAULT_MAXLOWDIM 8
#define DEFAULT_WINSIZE 25.0
#define DEFAULT_HOPSIZE 0.5
#define DEFAULT_DOWNSAMPLE 1
#define DEFAULT_MAXWINSIZE 500.0
#define DEFAULT_NORMALIZE 0
#define STR_(v) #v
#define STR(v) STR_(v)

// The sizes the analysis state is allocated for.
struct StateConfig {
    size_t highDim;
    size_t maxLowDim;
    double maxWinMs;
};

// C++ state lives behind a pointer because Max allocates objects with object_alloc (no constructors).
struct CccRpcState {
    const StateConfig config;

    double sampleRate = 0;
    size_t maxWindowSamples = 0;
    // maxLowDim x highDim, row-major. A projection into l dimensions uses the
    // first l rows, so @lowdim can change without touching the matrix. The
    // histogram rescales each dimension to its own range, so the matrix's
    // 1/sqrt(maxLowDim) scaling does not change the result.
    std::vector<double> matrix;
    std::vector<double> ringStorage;
    cccrt::RingBuffer<double> ring;
    std::vector<double> window;
    std::vector<double> projScratch;
    std::vector<uint64_t> cellScratch;
    size_t hopCounter = 0;
    // downsampling accumulator
    double accum = 0;
    size_t accumCount = 0;
    // scratch for noise calibration (@normalize 2), used on the main thread only
    std::vector<double> calibWindow;
    std::vector<double> calibProj;
    std::vector<uint64_t> calibCells;

    explicit CccRpcState(const StateConfig& c) : config(c) {
        matrix.resize(config.maxLowDim * config.highDim);
        cccrt::rpc::makeProjectionMatrix(matrix.data(), config.maxLowDim, config.highDim, 42);
    }

    bool matches(const StateConfig& c) const {
        return c.highDim == config.highDim && c.maxLowDim == config.maxLowDim && c.maxWinMs == config.maxWinMs;
    }

    // Called on the main thread (dsp64 or a rebuild), so allocation is fine here.
    void prepare(double sr) {
        if (sr == sampleRate && maxWindowSamples > 0) return;
        sampleRate = sr;
        maxWindowSamples = std::max<size_t>(1, static_cast<size_t>(config.maxWinMs / 1000.0 * sr));
        ringStorage.assign(maxWindowSamples + 1, 0.0);
        ring.init(ringStorage.data(), ringStorage.size());
        window.assign(maxWindowSamples, 0.0);
        // scratch for the worst case: a projection hop of one sample over the largest window
        const size_t maxHops = std::max<size_t>(1, cccrt::rpc::numHops(maxWindowSamples, config.highDim, 1));
        projScratch.assign(config.maxLowDim * maxHops, 0.0);
        cellScratch.assign(maxHops, 0);
        calibWindow.assign(maxWindowSamples, 0.0);
        calibProj.assign(config.maxLowDim * maxHops, 0.0);
        calibCells.assign(maxHops, 0);
        hopCounter = 0;
        accum = 0;
        accumCount = 0;
    }
};

typedef struct _cccrpc {
    t_pxobject ob;
    // The state the perform routine uses. Only the audio thread swaps it while
    // DSP runs (taking `incoming`, parking the old one in `outgoing`); only the
    // main thread frees states, so main-thread readers never see one freed.
    std::atomic<CccRpcState*> state;
    std::atomic<CccRpcState*> incoming; // rebuilt on the main thread, waiting for the audio thread
    std::atomic<CccRpcState*> outgoing; // replaced by the audio thread, waiting to be freed
    double output;   // latest value after normalisation, written by the perform routine
    double noiseRef; // reference for @normalize 2
    void* out_float;
    void* clock;   // defers float output from the perform routine to the scheduler
    void* serviceQelem; // main-thread upkeep after a change: free, rebuild, recalibrate
    // attributes
    long highdim;
    long lowdim;
    long res;
    double rpchop;
    long maxlowdim;
    double winsize; // ms
    double hopsize; // fraction of winsize
    long downsample;
    double maxwinsize; // ms
    long normalize;
} t_cccrpc;

static t_class* cccrpc_class = nullptr;

// the sizes the attributes ask for; @lowdim above @maxlowdim grows the allocation
static StateConfig cccrpc_config(const t_cccrpc* x) {
    StateConfig c;
    c.highDim = static_cast<size_t>(std::max<long>(1, x->highdim));
    c.maxLowDim = static_cast<size_t>(std::max<long>(1, std::max(x->maxlowdim, x->lowdim)));
    c.maxWinMs = x->maxwinsize > 0 ? x->maxwinsize : DEFAULT_MAXWINSIZE;
    return c;
}

// the analysis parameters derived from the attributes, shared by the perform routine and calibration
struct AnalysisParams {
    size_t ds;
    size_t windowSamples;
    size_t hopSamples;
    size_t resolution;
    size_t lowDim;
    double rpcHop;
};

static AnalysisParams cccrpc_params(const t_cccrpc* x, const CccRpcState& st) {
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
    // a @lowdim beyond this state's allocation applies once the rebuilt state arrives
    p.lowDim = std::min(static_cast<size_t>(std::max<long>(1, x->lowdim)), st.config.maxLowDim);
    p.rpcHop = x->rpchop;
    return p;
}

void* cccrpc_new(t_symbol* s, long argc, t_atom* argv);
void cccrpc_service(t_cccrpc* x);
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

    CLASS_ATTR_LONG(c, "highdim", 0, t_cccrpc, highdim);
    CLASS_ATTR_FILTER_MIN(c, "highdim", 1);
    CLASS_ATTR_LABEL(c, "highdim", 0, "Projection Window (h, samples)");
    CLASS_ATTR_CATEGORY(c, "highdim", 0, "Projection");
    CLASS_ATTR_ORDER(c, "highdim", 0, "1");
    CLASS_ATTR_DEFAULT_SAVE(c, "highdim", 0, STR(DEFAULT_HIGHDIM));

    CLASS_ATTR_LONG(c, "lowdim", 0, t_cccrpc, lowdim);
    CLASS_ATTR_FILTER_MIN(c, "lowdim", 1);
    CLASS_ATTR_LABEL(c, "lowdim", 0, "Projection Dimensions (l)");
    CLASS_ATTR_CATEGORY(c, "lowdim", 0, "Projection");
    CLASS_ATTR_ORDER(c, "lowdim", 0, "2");
    CLASS_ATTR_DEFAULT_SAVE(c, "lowdim", 0, STR(DEFAULT_LOWDIM));

    CLASS_ATTR_LONG(c, "res", 0, t_cccrpc, res);
    CLASS_ATTR_FILTER_MIN(c, "res", 1);
    CLASS_ATTR_LABEL(c, "res", 0, "Histogram Resolution");
    CLASS_ATTR_CATEGORY(c, "res", 0, "Projection");
    CLASS_ATTR_ORDER(c, "res", 0, "3");
    CLASS_ATTR_DEFAULT_SAVE(c, "res", 0, STR(DEFAULT_RES));

    CLASS_ATTR_DOUBLE(c, "rpchop", 0, t_cccrpc, rpchop);
    CLASS_ATTR_FILTER_CLIP(c, "rpchop", 0.0, 1.0);
    CLASS_ATTR_LABEL(c, "rpchop", 0, "Projection Hop (fraction of highdim)");
    CLASS_ATTR_CATEGORY(c, "rpchop", 0, "Projection");
    CLASS_ATTR_ORDER(c, "rpchop", 0, "4");
    CLASS_ATTR_DEFAULT_SAVE(c, "rpchop", 0, STR(DEFAULT_RPCHOP));

    CLASS_ATTR_LONG(c, "maxlowdim", 0, t_cccrpc, maxlowdim);
    CLASS_ATTR_FILTER_MIN(c, "maxlowdim", 1);
    CLASS_ATTR_LABEL(c, "maxlowdim", 0, "Preallocated Dimensions");
    CLASS_ATTR_CATEGORY(c, "maxlowdim", 0, "Projection");
    CLASS_ATTR_ORDER(c, "maxlowdim", 0, "5");
    CLASS_ATTR_DEFAULT_SAVE(c, "maxlowdim", 0, STR(DEFAULT_MAXLOWDIM));

    CLASS_ATTR_DOUBLE(c, "winsize", 0, t_cccrpc, winsize);
    CLASS_ATTR_FILTER_MIN(c, "winsize", 0.0);
    CLASS_ATTR_LABEL(c, "winsize", 0, "Window Size (ms)");
    CLASS_ATTR_CATEGORY(c, "winsize", 0, "Analysis");
    CLASS_ATTR_ORDER(c, "winsize", 0, "1");
    CLASS_ATTR_DEFAULT_SAVE(c, "winsize", 0, STR(DEFAULT_WINSIZE));

    CLASS_ATTR_DOUBLE(c, "hopsize", 0, t_cccrpc, hopsize);
    CLASS_ATTR_FILTER_MIN(c, "hopsize", 0.0);
    CLASS_ATTR_LABEL(c, "hopsize", 0, "Hop Size (fraction of window)");
    CLASS_ATTR_CATEGORY(c, "hopsize", 0, "Analysis");
    CLASS_ATTR_ORDER(c, "hopsize", 0, "2");
    CLASS_ATTR_DEFAULT_SAVE(c, "hopsize", 0, STR(DEFAULT_HOPSIZE));

    CLASS_ATTR_LONG(c, "downsample", 0, t_cccrpc, downsample);
    CLASS_ATTR_FILTER_MIN(c, "downsample", 1);
    CLASS_ATTR_LABEL(c, "downsample", 0, "Downsample Factor");
    CLASS_ATTR_CATEGORY(c, "downsample", 0, "Analysis");
    CLASS_ATTR_ORDER(c, "downsample", 0, "3");
    CLASS_ATTR_DEFAULT_SAVE(c, "downsample", 0, STR(DEFAULT_DOWNSAMPLE));

    CLASS_ATTR_DOUBLE(c, "maxwinsize", 0, t_cccrpc, maxwinsize);
    CLASS_ATTR_FILTER_MIN(c, "maxwinsize", 1.0);
    CLASS_ATTR_LABEL(c, "maxwinsize", 0, "Maximum Window Size (ms)");
    CLASS_ATTR_CATEGORY(c, "maxwinsize", 0, "Analysis");
    CLASS_ATTR_ORDER(c, "maxwinsize", 0, "4");
    CLASS_ATTR_DEFAULT_SAVE(c, "maxwinsize", 0, STR(DEFAULT_MAXWINSIZE));

    CLASS_ATTR_LONG(c, "normalize", 0, t_cccrpc, normalize);
    CLASS_ATTR_FILTER_CLIP(c, "normalize", 0, 2);
    CLASS_ATTR_ENUMINDEX3(c, "normalize", 0, "Raw", "Maximum", "White Noise");
    CLASS_ATTR_LABEL(c, "normalize", 0, "Normalize Output");
    CLASS_ATTR_CATEGORY(c, "normalize", 0, "Output");
    CLASS_ATTR_ORDER(c, "normalize", 0, "1");
    CLASS_ATTR_DEFAULT_SAVE(c, "normalize", 0, STR(DEFAULT_NORMALIZE));

    // every analysis parameter change re-runs the noise calibration, and a size change rebuilds the state
    for (const char* name : {"highdim", "lowdim", "res", "rpchop", "maxlowdim", "winsize", "hopsize", "downsample",
                             "maxwinsize"}) {
        CLASS_ATTR_ACCESSORS(c, name, NULL, cccrpc_attr_set);
    }

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    cccrpc_class = c;
}

void* cccrpc_new(t_symbol* s, long argc, t_atom* argv) {
    t_cccrpc* x = (t_cccrpc*)object_alloc(cccrpc_class);
    if (!x) return nullptr;
    new (&x->state) std::atomic<CccRpcState*>(nullptr);
    new (&x->incoming) std::atomic<CccRpcState*>(nullptr);
    new (&x->outgoing) std::atomic<CccRpcState*>(nullptr);
    x->output = 0;
    x->noiseRef = 1;

    x->highdim = DEFAULT_HIGHDIM;
    x->lowdim = DEFAULT_LOWDIM;
    x->res = DEFAULT_RES;
    x->rpchop = DEFAULT_RPCHOP;
    x->maxlowdim = DEFAULT_MAXLOWDIM;
    x->winsize = DEFAULT_WINSIZE;
    x->hopsize = DEFAULT_HOPSIZE;
    x->downsample = DEFAULT_DOWNSAMPLE;
    x->maxwinsize = DEFAULT_MAXWINSIZE;
    x->normalize = DEFAULT_NORMALIZE;

    // older patches: cccrpc~ [highDim] [lowDim] [maxWinSize] [maxLowDim]
    const long nPositional = attr_args_offset((short)argc, argv);
    if (nPositional > 0) x->highdim = std::max<long>(1, atom_getlong(argv));
    if (nPositional > 1) x->lowdim = std::max<long>(1, atom_getlong(argv + 1));
    if (nPositional > 2 && atom_getfloat(argv + 2) > 0.0) x->maxwinsize = atom_getfloat(argv + 2);
    if (nPositional > 3) x->maxlowdim = std::max<long>(1, atom_getlong(argv + 3));

    dsp_setup((t_pxobject*)x, 1);
    // outlets are created right to left
    x->out_float = floatout((t_object*)x);
    outlet_new((t_object*)x, "signal");
    x->clock = clock_new(x, (method)cccrpc_tick);
    x->serviceQelem = qelem_new(x, (method)cccrpc_service);

    // typed attributes; the setters only store values until the state exists
    attr_args_process(x, (short)argc, argv);
    x->state.store(new CccRpcState(cccrpc_config(x)));
    return x;
}

void cccrpc_free(t_cccrpc* x) {
    dsp_free((t_pxobject*)x); // the perform routine no longer runs after this
    object_free(x->clock);
    qelem_free(x->serviceQelem);
    delete x->state.exchange(nullptr);
    delete x->incoming.exchange(nullptr);
    delete x->outgoing.exchange(nullptr);
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
    outlet_float(x->out_float, x->output);
}

// shared setter for the analysis attributes: store the value, then let the service qelem catch up
t_max_err cccrpc_attr_set(t_cccrpc* x, t_object* attr, long argc, t_atom* argv) {
    if (argc < 1 || !argv) return MAX_ERR_NONE;
    t_symbol* name = (t_symbol*)object_method(attr, gensym("getname"));
    if (name == gensym("highdim")) x->highdim = std::max<long>(1, atom_getlong(argv));
    else if (name == gensym("lowdim")) x->lowdim = std::max<long>(1, atom_getlong(argv));
    else if (name == gensym("res")) x->res = std::max<long>(1, atom_getlong(argv));
    else if (name == gensym("rpchop")) x->rpchop = std::min(1.0, std::max(0.0, atom_getfloat(argv)));
    else if (name == gensym("maxlowdim")) x->maxlowdim = std::max<long>(1, atom_getlong(argv));
    else if (name == gensym("winsize")) x->winsize = std::max(0.0, atom_getfloat(argv));
    else if (name == gensym("hopsize")) x->hopsize = std::max(0.0, atom_getfloat(argv));
    else if (name == gensym("downsample")) x->downsample = std::max<long>(1, atom_getlong(argv));
    else if (name == gensym("maxwinsize")) x->maxwinsize = std::max(1.0, atom_getfloat(argv));
    if (x->serviceQelem) qelem_set(x->serviceQelem);
    return MAX_ERR_NONE;
}

// Main-thread upkeep, run by the qelem after an attribute change or a state
// swap: free a state the audio thread has retired, rebuild if the sizes the
// attributes ask for have changed, and recalibrate the noise reference.
// Several attribute changes in a row (e.g. loading a patch) cost one rebuild.
void cccrpc_service(t_cccrpc* x) {
    delete x->outgoing.exchange(nullptr);
    CccRpcState* current = x->state.load();
    if (!current) return;

    const StateConfig want = cccrpc_config(x);
    CccRpcState* pending = x->incoming.load();
    if (!(pending ? pending : current)->matches(want)) {
        CccRpcState* next = new CccRpcState(want);
        if (current->sampleRate > 0) next->prepare(current->sampleRate);
        if (sys_getdspobjdspstate((t_object*)x)) {
            // the perform routine takes it at its next vector, and the qelem runs again to tidy up
            delete x->incoming.exchange(next);
        } else {
            // no perform routine running: swap in place
            delete x->incoming.exchange(nullptr);
            x->state.store(next);
            delete current;
        }
    }
    cccrpc_calibrate(x);
}

// Reference value for @normalize 2: mean RPC of white noise windows with the
// current parameters. Runs on the main thread (qelem / dsp64) with its own scratch.
void cccrpc_calibrate(t_cccrpc* x) {
    CccRpcState& st = *x->state.load();
    if (st.sampleRate <= 0 || st.calibWindow.empty()) return; // dsp64 will calibrate once it knows the sample rate
    const AnalysisParams p = cccrpc_params(x, st);
    cccrt::rpc::Pcg32 rng(1234); // fixed seed: the same parameters always give the same reference
    const int nWindows = 8;
    double sum = 0;
    for (int w = 0; w < nWindows; ++w) {
        for (size_t i = 0; i < p.windowSamples; ++i) st.calibWindow[i] = rng.uniform01<double>() * 2.0 - 1.0;
        sum += cccrt::rpc::calc(st.matrix.data(), p.lowDim, st.config.highDim, st.calibWindow.data(),
                                p.windowSamples, p.resolution, p.rpcHop, st.calibProj.data(), st.calibCells.data());
    }
    x->noiseRef = std::max(1.0, sum / nWindows);
}

void cccrpc_dsp64(t_cccrpc* x, t_object* dsp64, short* count, double samplerate, long maxvectorsize, long flags) {
    x->state.load()->prepare(samplerate);
    // a waiting state was prepared for the rate current when it was built; a
    // new rate means audio restarted, so nothing is using it yet
    if (CccRpcState* pending = x->incoming.load()) pending->prepare(samplerate);
    cccrpc_calibrate(x);
    object_method(dsp64, gensym("dsp_add64"), x, cccrpc_perform64, 0, NULL);
}

void cccrpc_perform64(t_cccrpc* x, t_object* dsp64, double** ins, long numins, double** outs, long numouts,
                      long sampleframes, long flags, void* userparam) {
    // pick up a rebuilt state, once the previous swap has been tidied up
    if (!x->outgoing.load(std::memory_order_acquire)) {
        if (CccRpcState* next = x->incoming.exchange(nullptr, std::memory_order_acq_rel)) {
            x->outgoing.store(x->state.exchange(next, std::memory_order_acq_rel), std::memory_order_release);
            qelem_set(x->serviceQelem);
        }
    }
    CccRpcState& st = *x->state.load(std::memory_order_acquire);
    const double* in = ins[0];
    double* out = outs[0];

    // read the modulatable parameters once per block, as the UGen does
    const AnalysisParams p = cccrpc_params(x, st);
    const size_t ds = p.ds;
    const size_t windowSamples = p.windowSamples;
    const size_t hopSamples = p.hopSamples;
    const size_t resolution = p.resolution;
    const double rpcHop = p.rpcHop;
    const size_t lowDim = p.lowDim;
    const long normalize = x->normalize;
    double output = x->output;
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
                const double rpc = cccrt::rpc::calc(st.matrix.data(), lowDim, st.config.highDim, st.window.data(),
                                                    windowSamples, resolution, rpcHop, st.projScratch.data(),
                                                    st.cellScratch.data());
                if (normalize == 1) {
                    const double mx =
                        cccrt::rpc::maxOccupiedCells(windowSamples, st.config.highDim, rpcHop, resolution, lowDim);
                    output = mx > 0 ? rpc / mx : 0;
                } else if (normalize == 2) {
                    output = rpc / x->noiseRef;
                } else {
                    output = rpc;
                }
                newValue = true;
            }
        }
        out[i] = output;
    }
    x->output = output;
    if (newValue) clock_delay(x->clock, 0);
}
