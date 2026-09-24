// cccrpc - Random Projection Complexity of a frame of data (control rate)
// Chris Kiefer (c.kiefer@sussex.ac.uk)
//
// The control-rate companion to cccrpc~. Where cccrpc~ analyses a sliding
// window of an audio signal, this analyses one complete frame handed to it as
// a list or read from a buffer~ - an FFT magnitude spectrum, a sensor frame, a
// table of values. Each frame produces one value.
//
// Fed a spectrum, RPC slides its projection window along the FREQUENCY axis:
// @highdim is a number of adjacent bins and @rpchop is a bin hop. A spectrum
// with regular structure (a harmonic series) projects repeatedly onto the same
// histogram cells and scores low; a dense or noisy spectrum scores high.
//
//   cccrpc [@attribute value ...]
//   inlet
//     list       : analyse these values as one frame, output the result
//     bang       : read the buffer~ named by @buffer and analyse that
//     buffer <n> : set the buffer~ name (same as @buffer)
//   outlets
//     0 (float)  : the complexity value
//     1 (int)    : how many values were analysed after @skip / @bins
//   attributes
//     @buffer    : buffer~ to read on bang
//     @skip      : ignore this many values at the start of the frame (default 0;
//                  set 1 for an FFT spectrum to drop the DC bin)
//     @bins      : use at most this many values after @skip (0 = all, the default;
//                  e.g. 256 of 512 bins to ignore the top octave)
//     @logmag    : 1 = convert values to dB before analysis (default 0). Worth
//                  turning on for FFT magnitudes: a few large low-frequency bins
//                  would otherwise squash the rest into one histogram cell.
//     @dbfloor   : floor for @logmag, in dB (default -120)
//     @maxframe  : longest frame accepted, in values (default 4096)
//     @highdim   : projection window length in values, h (default 16)
//     @lowdim    : projection dimensions, l (default 2)
//     @res       : histogram resolution per dimension, beta (default 10)
//     @rpchop    : projection hop as a fraction of highdim (default 0.5)
//     @maxlowdim : @lowdim values up to this are preallocated (default 8); a
//                  larger @lowdim still works, but reallocates
//     @normalize : 0 = raw cell count (default)
//                  1 = divided by the maximum possible count, min(hops, res^l)
//                  2 = divided by the mean count for a random frame of the same
//                      length (so noise reads about 1 whatever the settings)
//
//   Older patches may give cccrpc [highDim] [lowDim] [maxFrame] [maxLowDim]
//   as arguments; these set the matching attributes, and a typed or saved
//   attribute takes precedence.
//
// Changing @highdim, @maxframe, or @lowdim beyond the allocated size rebuilds
// the analysis state. A critical region keeps the rebuild (main thread) from
// overlapping an analysis (main or scheduler thread, depending on overdrive).

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "ext.h"
#include "ext_obex.h"
#include "ext_buffer.h"
#include "ext_critical.h"

#include "rpc.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

// Attribute defaults. The STR() of each feeds CLASS_ATTR_DEFAULT_SAVE, so the
// inspector and cccrpc_new agree (Max does not apply attribute defaults to
// non-UI objects; cccrpc_new writes them into the struct).
#define DEFAULT_SKIP 0
#define DEFAULT_BINS 0
#define DEFAULT_LOGMAG 0
#define DEFAULT_DBFLOOR -120.0
#define DEFAULT_MAXFRAME 4096
#define DEFAULT_HIGHDIM 16
#define DEFAULT_LOWDIM 2
#define DEFAULT_RES 10
#define DEFAULT_RPCHOP 0.5
#define DEFAULT_MAXLOWDIM 8
#define DEFAULT_NORMALIZE 0
#define STR_(v) #v
#define STR(v) STR_(v)

// The sizes the analysis state is allocated for.
struct StateConfig {
    size_t highDim;
    size_t maxLowDim;
    size_t maxFrame;
};

// C++ state lives behind a pointer because Max allocates objects with object_alloc (no constructors).
struct CccRpcState {
    const StateConfig config;

    // maxLowDim x highDim, row-major. A projection into l dimensions uses the
    // first l rows, so @lowdim can change without rebuilding the matrix. The
    // histogram rescales each dimension to its own range, so the matrix's
    // 1/sqrt(maxLowDim) scaling does not change the result.
    std::vector<double> matrix;
    std::vector<double> raw;     // frame as received
    std::vector<double> frame;   // after @skip / @bins / @logmag
    std::vector<double> proj;
    std::vector<uint64_t> cells;
    // noise calibration for @normalize 2, cached against the settings it was made with
    std::vector<double> calibFrame;
    std::vector<double> calibProj;
    std::vector<uint64_t> calibCells;
    double noiseRef = 1;
    size_t refN = 0, refLowDim = 0, refRes = 0;
    double refHop = -1;

    explicit CccRpcState(const StateConfig& c) : config(c) {
        matrix.resize(config.maxLowDim * config.highDim);
        cccrt::rpc::makeProjectionMatrix(matrix.data(), config.maxLowDim, config.highDim, 42);
        raw.assign(config.maxFrame, 0.0);
        frame.assign(config.maxFrame, 0.0);
        calibFrame.assign(config.maxFrame, 0.0);
        // worst case: a projection hop of one value over the longest frame
        const size_t maxHops = std::max<size_t>(1, cccrt::rpc::numHops(config.maxFrame, config.highDim, 1));
        proj.assign(config.maxLowDim * maxHops, 0.0);
        cells.assign(maxHops, 0);
        calibProj.assign(config.maxLowDim * maxHops, 0.0);
        calibCells.assign(maxHops, 0);
    }

    bool matches(const StateConfig& c) const {
        return c.highDim == config.highDim && c.maxLowDim == config.maxLowDim && c.maxFrame == config.maxFrame;
    }
};

typedef struct _cccrpc {
    t_object ob;
    CccRpcState* state; // swapped under `lock`
    t_critical lock;
    t_buffer_ref* bufref;
    t_symbol* buffername;
    void* out_count;   // right outlet
    void* out_value;   // left outlet
    // attributes
    long skip;
    long bins;
    long logmag;
    double dbfloor;
    long maxframe;
    long highdim;
    long lowdim;
    long res;
    double rpchop;
    long maxlowdim;
    long normalize;
} t_cccrpc;

static t_class* cccrpc_class = nullptr;

void* cccrpc_new(t_symbol* s, long argc, t_atom* argv);
void cccrpc_free(t_cccrpc* x);
void cccrpc_assist(t_cccrpc* x, void* b, long m, long a, char* s);
void cccrpc_list(t_cccrpc* x, t_symbol* s, long argc, t_atom* argv);
void cccrpc_bang(t_cccrpc* x);
void cccrpc_buffer(t_cccrpc* x, t_symbol* name);
t_max_err cccrpc_notify(t_cccrpc* x, t_symbol* s, t_symbol* msg, void* sender, void* data);
t_max_err cccrpc_set_buffer(t_cccrpc* x, t_object* attr, long argc, t_atom* argv);
t_max_err cccrpc_set_size(t_cccrpc* x, t_object* attr, long argc, t_atom* argv);

void C74_EXPORT ext_main(void* r) {
    t_class* c = class_new("cccrpc", (method)cccrpc_new, (method)cccrpc_free, (long)sizeof(t_cccrpc), 0L, A_GIMME, 0);

    class_addmethod(c, (method)cccrpc_list, "list", A_GIMME, 0);
    class_addmethod(c, (method)cccrpc_bang, "bang", 0);
    class_addmethod(c, (method)cccrpc_buffer, "buffer", A_SYM, 0);
    class_addmethod(c, (method)cccrpc_notify, "notify", A_CANT, 0);
    class_addmethod(c, (method)cccrpc_assist, "assist", A_CANT, 0);

    CLASS_ATTR_SYM(c, "buffer", 0, t_cccrpc, buffername);
    CLASS_ATTR_ACCESSORS(c, "buffer", NULL, cccrpc_set_buffer);
    CLASS_ATTR_LABEL(c, "buffer", 0, "buffer~ to analyse on bang");
    CLASS_ATTR_CATEGORY(c, "buffer", 0, "Input");
    CLASS_ATTR_ORDER(c, "buffer", 0, "1");
    CLASS_ATTR_SAVE(c, "buffer", 0);

    CLASS_ATTR_LONG(c, "skip", 0, t_cccrpc, skip);
    CLASS_ATTR_FILTER_MIN(c, "skip", 0);
    CLASS_ATTR_LABEL(c, "skip", 0, "Values To Skip");
    CLASS_ATTR_CATEGORY(c, "skip", 0, "Input");
    CLASS_ATTR_ORDER(c, "skip", 0, "2");
    CLASS_ATTR_DEFAULT_SAVE(c, "skip", 0, STR(DEFAULT_SKIP));

    CLASS_ATTR_LONG(c, "bins", 0, t_cccrpc, bins);
    CLASS_ATTR_FILTER_MIN(c, "bins", 0);
    CLASS_ATTR_LABEL(c, "bins", 0, "Values To Use (0 = all)");
    CLASS_ATTR_CATEGORY(c, "bins", 0, "Input");
    CLASS_ATTR_ORDER(c, "bins", 0, "3");
    CLASS_ATTR_DEFAULT_SAVE(c, "bins", 0, STR(DEFAULT_BINS));

    CLASS_ATTR_LONG(c, "logmag", 0, t_cccrpc, logmag);
    CLASS_ATTR_STYLE_LABEL(c, "logmag", 0, "onoff", "Convert To dB");
    CLASS_ATTR_CATEGORY(c, "logmag", 0, "Input");
    CLASS_ATTR_ORDER(c, "logmag", 0, "4");
    CLASS_ATTR_DEFAULT_SAVE(c, "logmag", 0, STR(DEFAULT_LOGMAG));

    CLASS_ATTR_DOUBLE(c, "dbfloor", 0, t_cccrpc, dbfloor);
    CLASS_ATTR_LABEL(c, "dbfloor", 0, "dB Floor For logmag");
    CLASS_ATTR_CATEGORY(c, "dbfloor", 0, "Input");
    CLASS_ATTR_ORDER(c, "dbfloor", 0, "5");
    CLASS_ATTR_DEFAULT_SAVE(c, "dbfloor", 0, STR(DEFAULT_DBFLOOR));

    CLASS_ATTR_LONG(c, "maxframe", 0, t_cccrpc, maxframe);
    CLASS_ATTR_ACCESSORS(c, "maxframe", NULL, cccrpc_set_size);
    CLASS_ATTR_LABEL(c, "maxframe", 0, "Maximum Frame Length");
    CLASS_ATTR_CATEGORY(c, "maxframe", 0, "Input");
    CLASS_ATTR_ORDER(c, "maxframe", 0, "6");
    CLASS_ATTR_DEFAULT_SAVE(c, "maxframe", 0, STR(DEFAULT_MAXFRAME));

    CLASS_ATTR_LONG(c, "highdim", 0, t_cccrpc, highdim);
    CLASS_ATTR_ACCESSORS(c, "highdim", NULL, cccrpc_set_size);
    CLASS_ATTR_LABEL(c, "highdim", 0, "Projection Window (h, values)");
    CLASS_ATTR_CATEGORY(c, "highdim", 0, "Projection");
    CLASS_ATTR_ORDER(c, "highdim", 0, "1");
    CLASS_ATTR_DEFAULT_SAVE(c, "highdim", 0, STR(DEFAULT_HIGHDIM));

    CLASS_ATTR_LONG(c, "lowdim", 0, t_cccrpc, lowdim);
    CLASS_ATTR_ACCESSORS(c, "lowdim", NULL, cccrpc_set_size);
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
    CLASS_ATTR_ACCESSORS(c, "maxlowdim", NULL, cccrpc_set_size);
    CLASS_ATTR_LABEL(c, "maxlowdim", 0, "Preallocated Dimensions");
    CLASS_ATTR_CATEGORY(c, "maxlowdim", 0, "Projection");
    CLASS_ATTR_ORDER(c, "maxlowdim", 0, "5");
    CLASS_ATTR_DEFAULT_SAVE(c, "maxlowdim", 0, STR(DEFAULT_MAXLOWDIM));

    CLASS_ATTR_LONG(c, "normalize", 0, t_cccrpc, normalize);
    CLASS_ATTR_FILTER_CLIP(c, "normalize", 0, 2);
    CLASS_ATTR_ENUMINDEX3(c, "normalize", 0, "Raw", "Maximum", "Random Frame");
    CLASS_ATTR_LABEL(c, "normalize", 0, "Normalize Output");
    CLASS_ATTR_CATEGORY(c, "normalize", 0, "Output");
    CLASS_ATTR_ORDER(c, "normalize", 0, "1");
    CLASS_ATTR_DEFAULT_SAVE(c, "normalize", 0, STR(DEFAULT_NORMALIZE));

    class_register(CLASS_BOX, c);
    cccrpc_class = c;
}

// the sizes the attributes ask for; @lowdim above @maxlowdim grows the allocation
static StateConfig cccrpc_config(const t_cccrpc* x) {
    StateConfig c;
    c.highDim = static_cast<size_t>(std::max<long>(1, x->highdim));
    c.maxLowDim = static_cast<size_t>(std::max<long>(1, std::max(x->maxlowdim, x->lowdim)));
    c.maxFrame = static_cast<size_t>(std::max(x->highdim, x->maxframe));
    return c;
}

// Reallocate if the size attributes have changed. Main thread only.
static void cccrpc_rebuild(t_cccrpc* x) {
    if (!x->state) return; // still in cccrpc_new, which builds it after the attributes
    const StateConfig want = cccrpc_config(x);
    if (x->state->matches(want)) return;
    CccRpcState* next = new CccRpcState(want);
    critical_enter(x->lock);
    CccRpcState* old = x->state;
    x->state = next;
    critical_exit(x->lock);
    delete old;
}

void* cccrpc_new(t_symbol* s, long argc, t_atom* argv) {
    t_cccrpc* x = (t_cccrpc*)object_alloc(cccrpc_class);
    if (!x) return nullptr;

    x->state = nullptr;
    x->bufref = nullptr;
    x->buffername = gensym("");
    x->skip = DEFAULT_SKIP;
    x->bins = DEFAULT_BINS;
    x->logmag = DEFAULT_LOGMAG;
    x->dbfloor = DEFAULT_DBFLOOR;
    x->maxframe = DEFAULT_MAXFRAME;
    x->highdim = DEFAULT_HIGHDIM;
    x->lowdim = DEFAULT_LOWDIM;
    x->res = DEFAULT_RES;
    x->rpchop = DEFAULT_RPCHOP;
    x->maxlowdim = DEFAULT_MAXLOWDIM;
    x->normalize = DEFAULT_NORMALIZE;

    // older patches: cccrpc [highDim] [lowDim] [maxFrame] [maxLowDim]
    const long nPositional = attr_args_offset((short)argc, argv);
    if (nPositional > 0) x->highdim = std::max<long>(1, atom_getlong(argv));
    if (nPositional > 1) x->lowdim = std::max<long>(1, atom_getlong(argv + 1));
    if (nPositional > 2) x->maxframe = std::max<long>(1, atom_getlong(argv + 2));
    if (nPositional > 3) x->maxlowdim = std::max<long>(1, atom_getlong(argv + 3));

    critical_new(&x->lock);

    // outlets are created right to left
    x->out_count = intout((t_object*)x);
    x->out_value = floatout((t_object*)x);

    // typed attributes; the size setters only store values until the state exists
    attr_args_process(x, (short)argc, argv);
    x->state = new CccRpcState(cccrpc_config(x));
    return x;
}

void cccrpc_free(t_cccrpc* x) {
    object_free(x->bufref);
    delete x->state;
    x->state = nullptr;
    critical_free(x->lock);
}

void cccrpc_assist(t_cccrpc* x, void* b, long m, long a, char* s) {
    if (m == ASSIST_INLET) {
        snprintf(s, 256, "(list) Frame to analyse, (bang) analyse buffer~");
    } else if (a == 0) {
        snprintf(s, 256, "(float) Random projection complexity");
    } else {
        snprintf(s, 256, "(int) Number of values analysed");
    }
}

void cccrpc_buffer(t_cccrpc* x, t_symbol* name) {
    x->buffername = name;
    if (x->bufref)
        buffer_ref_set(x->bufref, name);
    else
        x->bufref = buffer_ref_new((t_object*)x, name);
}

t_max_err cccrpc_set_buffer(t_cccrpc* x, t_object* attr, long argc, t_atom* argv) {
    if (argc > 0 && argv) cccrpc_buffer(x, atom_getsym(argv));
    return MAX_ERR_NONE;
}

// setter for the attributes that decide how much the state allocates
t_max_err cccrpc_set_size(t_cccrpc* x, t_object* attr, long argc, t_atom* argv) {
    if (argc < 1 || !argv) return MAX_ERR_NONE;
    t_symbol* name = (t_symbol*)object_method(attr, gensym("getname"));
    const long v = std::max<long>(1, atom_getlong(argv));
    if (name == gensym("maxframe")) x->maxframe = v;
    else if (name == gensym("highdim")) x->highdim = v;
    else if (name == gensym("lowdim")) x->lowdim = v;
    else if (name == gensym("maxlowdim")) x->maxlowdim = v;
    cccrpc_rebuild(x);
    return MAX_ERR_NONE;
}

t_max_err cccrpc_notify(t_cccrpc* x, t_symbol* s, t_symbol* msg, void* sender, void* data) {
    return x->bufref ? buffer_ref_notify(x->bufref, s, msg, sender, data) : MAX_ERR_NONE;
}

// Apply @skip / @bins / @logmag to state.raw, writing state.frame. Returns the length.
static size_t cccrpc_prepare(t_cccrpc* x, CccRpcState& st, size_t n) {
    const size_t start = std::min(static_cast<size_t>(std::max<long>(0, x->skip)), n);
    size_t count = n - start;
    if (x->bins > 0) count = std::min(count, static_cast<size_t>(x->bins));
    count = std::min(count, st.config.maxFrame);

    if (x->logmag) {
        const double floorAmp = std::pow(10.0, x->dbfloor / 20.0);
        for (size_t i = 0; i < count; ++i) {
            st.frame[i] = 20.0 * std::log10(std::max(std::abs(st.raw[start + i]), floorAmp));
        }
    } else {
        std::copy(st.raw.begin() + start, st.raw.begin() + start + count, st.frame.begin());
    }
    return count;
}

// Mean RPC of random frames of this length with these settings, cached.
// Everything here runs on the main/scheduler thread, so the cost is not in an audio path.
static double cccrpc_noise_reference(CccRpcState& st, size_t n, size_t lowDim, size_t res, double hop) {
    if (st.refN == n && st.refLowDim == lowDim && st.refRes == res && st.refHop == hop) return st.noiseRef;

    cccrt::rpc::Pcg32 rng(1234);   // fixed seed: the same settings always give the same reference
    const int nFrames = 8;
    double sum = 0;
    for (int f = 0; f < nFrames; ++f) {
        for (size_t i = 0; i < n; ++i) st.calibFrame[i] = rng.uniform01<double>() * 2.0 - 1.0;
        sum += cccrt::rpc::calc(st.matrix.data(), lowDim, st.config.highDim, st.calibFrame.data(), n,
                                res, hop, st.calibProj.data(), st.calibCells.data());
    }
    st.noiseRef = std::max(1.0, sum / nFrames);
    st.refN = n; st.refLowDim = lowDim; st.refRes = res; st.refHop = hop;
    return st.noiseRef;
}

// Analyse state.frame[0..count). Called with the lock held.
static double cccrpc_analyse(t_cccrpc* x, CccRpcState& st, size_t count) {
    if (count < st.config.highDim) return 0.0;   // too short to project even once
    const size_t lowDim = std::min(static_cast<size_t>(std::max<long>(1, x->lowdim)), st.config.maxLowDim);
    const size_t res = static_cast<size_t>(std::max<long>(1, x->res));
    const double hop = x->rpchop;

    double value = cccrt::rpc::calc(st.matrix.data(), lowDim, st.config.highDim, st.frame.data(), count,
                                    res, hop, st.proj.data(), st.cells.data());
    if (x->normalize == 1) {
        const double mx = cccrt::rpc::maxOccupiedCells(count, st.config.highDim, hop, res, lowDim);
        value = mx > 0 ? value / mx : 0;
    } else if (x->normalize == 2) {
        value /= cccrpc_noise_reference(st, count, lowDim, res, hop);
    }
    return value;
}

// outlets fire after the lock is released, so downstream objects can set our attributes
static void cccrpc_output(t_cccrpc* x, size_t count, double value) {
    outlet_int(x->out_count, static_cast<t_atom_long>(count));
    outlet_float(x->out_value, value);
}

void cccrpc_list(t_cccrpc* x, t_symbol* s, long argc, t_atom* argv) {
    critical_enter(x->lock);
    CccRpcState& st = *x->state;
    const size_t maxFrame = st.config.maxFrame;
    const size_t n = std::min(static_cast<size_t>(std::max<long>(0, argc)), maxFrame);
    for (size_t i = 0; i < n; ++i) st.raw[i] = atom_getfloat(argv + i);
    const size_t count = cccrpc_prepare(x, st, n);
    const double value = cccrpc_analyse(x, st, count);
    critical_exit(x->lock);

    if (static_cast<size_t>(argc) > maxFrame) {
        object_warn((t_object*)x, "frame of %ld values truncated to @maxframe (%ld)", argc, (long)maxFrame);
    }
    cccrpc_output(x, count, value);
}

void cccrpc_bang(t_cccrpc* x) {
    if (!x->bufref) {
        object_error((t_object*)x, "no buffer~ set: send 'buffer <name>' or use @buffer");
        return;
    }
    t_buffer_obj* buffer = buffer_ref_getobject(x->bufref);
    if (!buffer) {
        object_error((t_object*)x, "no buffer~ named '%s'", x->buffername->s_name);
        return;
    }
    float* samples = buffer_locksamples(buffer);
    if (!samples) {
        object_error((t_object*)x, "buffer~ '%s' is busy", x->buffername->s_name);
        return;
    }
    const size_t channels = std::max<size_t>(1, static_cast<size_t>(buffer_getchannelcount(buffer)));
    const size_t frames = static_cast<size_t>(buffer_getframecount(buffer));

    critical_enter(x->lock);
    CccRpcState& st = *x->state;
    const size_t n = std::min(frames, st.config.maxFrame);
    for (size_t i = 0; i < n; ++i) st.raw[i] = samples[i * channels];   // channel 0
    buffer_unlocksamples(buffer);
    const size_t count = cccrpc_prepare(x, st, n);
    const double value = cccrpc_analyse(x, st, count);
    critical_exit(x->lock);

    cccrpc_output(x, count, value);
}
