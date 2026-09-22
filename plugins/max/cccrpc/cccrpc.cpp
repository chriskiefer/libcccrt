// cccrpc - Random Projection Complexity of a frame of data (control rate)
// Chris Kiefer (c.kiefer@sussex.ac.uk)
//
// The control-rate companion to cccrpc~. Where cccrpc~ analyses a sliding
// window of an audio signal, this analyses one complete frame handed to it as
// a list or read from a buffer~ - an FFT magnitude spectrum, a sensor frame, a
// table of values. Each frame produces one value.
//
// Fed a spectrum, RPC slides its projection window along the FREQUENCY axis:
// highDim is a number of adjacent bins and @rpchop is a bin hop. A spectrum
// with regular structure (a harmonic series) projects repeatedly onto the same
// histogram cells and scores low; a dense or noisy spectrum scores high.
//
//   cccrpc [highDim] [lowDim] [maxFrame] [maxLowDim]
//     highDim    : projection window length in values, h (default 16) - fixed at creation
//     lowDim     : initial projection dimensions, l (default 2)
//     maxFrame   : longest frame accepted, in values (default 4096)   - fixed at creation
//     maxLowDim  : upper limit for @lowdim (default 8, or lowDim if larger) - fixed at creation
//   inlet
//     list       : analyse these values as one frame, output the result
//     bang       : read the buffer~ named by @buffer and analyse that
//     buffer <n> : set the buffer~ name (same as @buffer)
//   outlets
//     0 (float)  : the complexity value
//     1 (int)    : how many values were analysed after @skip / @bins
//   attributes
//     @buffer    : buffer~ to read on bang
//     @lowdim    : projection dimensions, l (1 .. maxLowDim)
//     @res       : histogram resolution per dimension, beta (default 10)
//     @rpchop    : projection hop as a fraction of highDim (default 0.5)
//     @normalize : 0 = raw cell count (default)
//                  1 = divided by the maximum possible count, min(hops, res^l)
//                  2 = divided by the mean count for a random frame of the same
//                      length (so noise reads about 1 whatever the settings)
//     @skip      : ignore this many values at the start of the frame (default 0;
//                  set 1 for an FFT spectrum to drop the DC bin)
//     @bins      : use at most this many values after @skip (0 = all, the default;
//                  e.g. 256 of 512 bins to ignore the top octave)
//     @logmag    : 1 = convert values to dB before analysis (default 0). Worth
//                  turning on for FFT magnitudes: a few large low-frequency bins
//                  would otherwise squash the rest into one histogram cell.
//     @dbfloor   : floor for @logmag, in dB (default -120)

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "ext.h"
#include "ext_obex.h"
#include "ext_buffer.h"

#include "rpc.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

// C++ state lives behind a pointer because Max allocates objects with object_alloc (no constructors).
struct CccRpcState {
    const size_t highDim;
    const size_t maxLowDim;
    const size_t maxFrame;

    // maxLowDim x highDim, row-major. A projection into l dimensions uses the
    // first l rows, so @lowdim can change without rebuilding the matrix.
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

    CccRpcState(size_t hd, size_t maxLd, size_t maxN)
        : highDim(hd), maxLowDim(maxLd), maxFrame(maxN) {
        matrix.resize(maxLowDim * highDim);
        cccrt::rpc::makeProjectionMatrix(matrix.data(), maxLowDim, highDim, 42);
        raw.assign(maxFrame, 0.0);
        frame.assign(maxFrame, 0.0);
        calibFrame.assign(maxFrame, 0.0);
        // worst case: a projection hop of one value over the longest frame
        const size_t maxHops = std::max<size_t>(1, cccrt::rpc::numHops(maxFrame, highDim, 1));
        proj.assign(maxLowDim * maxHops, 0.0);
        cells.assign(maxHops, 0);
        calibProj.assign(maxLowDim * maxHops, 0.0);
        calibCells.assign(maxHops, 0);
    }
};

typedef struct _cccrpc {
    t_object ob;
    CccRpcState* state;
    t_buffer_ref* bufref;
    t_symbol* buffername;
    void* out_count;   // right outlet
    void* out_value;   // left outlet
    // attributes
    long lowdim;
    long res;
    double rpchop;
    long normalize;
    long skip;
    long bins;
    long logmag;
    double dbfloor;
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

    CLASS_ATTR_LONG(c, "lowdim", 0, t_cccrpc, lowdim);
    CLASS_ATTR_FILTER_MIN(c, "lowdim", 1);
    CLASS_ATTR_LABEL(c, "lowdim", 0, "Projection Dimensions (l)");

    CLASS_ATTR_LONG(c, "res", 0, t_cccrpc, res);
    CLASS_ATTR_FILTER_MIN(c, "res", 1);
    CLASS_ATTR_LABEL(c, "res", 0, "Histogram Resolution");

    CLASS_ATTR_DOUBLE(c, "rpchop", 0, t_cccrpc, rpchop);
    CLASS_ATTR_FILTER_CLIP(c, "rpchop", 0.0, 1.0);
    CLASS_ATTR_LABEL(c, "rpchop", 0, "Projection Hop (fraction of highDim)");

    CLASS_ATTR_LONG(c, "normalize", 0, t_cccrpc, normalize);
    CLASS_ATTR_FILTER_CLIP(c, "normalize", 0, 2);
    CLASS_ATTR_ENUMINDEX3(c, "normalize", 0, "Raw", "Maximum", "Random Frame");
    CLASS_ATTR_LABEL(c, "normalize", 0, "Normalize Output");

    CLASS_ATTR_LONG(c, "skip", 0, t_cccrpc, skip);
    CLASS_ATTR_FILTER_MIN(c, "skip", 0);
    CLASS_ATTR_LABEL(c, "skip", 0, "Values To Skip");

    CLASS_ATTR_LONG(c, "bins", 0, t_cccrpc, bins);
    CLASS_ATTR_FILTER_MIN(c, "bins", 0);
    CLASS_ATTR_LABEL(c, "bins", 0, "Values To Use (0 = all)");

    CLASS_ATTR_LONG(c, "logmag", 0, t_cccrpc, logmag);
    CLASS_ATTR_STYLE_LABEL(c, "logmag", 0, "onoff", "Convert To dB");

    CLASS_ATTR_DOUBLE(c, "dbfloor", 0, t_cccrpc, dbfloor);
    CLASS_ATTR_LABEL(c, "dbfloor", 0, "dB Floor For logmag");

    class_register(CLASS_BOX, c);
    cccrpc_class = c;
}

void* cccrpc_new(t_symbol* s, long argc, t_atom* argv) {
    t_cccrpc* x = (t_cccrpc*)object_alloc(cccrpc_class);
    if (!x) return nullptr;

    // positional args, then attributes
    const long nPositional = attr_args_offset((short)argc, argv);
    long highDim = 16;
    long lowDim = 2;
    long maxFrame = 4096;
    long maxLowDim = 8;
    if (nPositional > 0) highDim = atom_getlong(argv);
    if (nPositional > 1) lowDim = atom_getlong(argv + 1);
    if (nPositional > 2) maxFrame = atom_getlong(argv + 2);
    if (nPositional > 3) maxLowDim = atom_getlong(argv + 3);
    highDim = std::max<long>(1, highDim);
    lowDim = std::max<long>(1, lowDim);
    maxFrame = std::max(highDim, maxFrame);
    maxLowDim = std::max(std::max<long>(1, maxLowDim), lowDim);

    x->bufref = nullptr;
    x->buffername = gensym("");
    x->lowdim = lowDim;
    x->res = 10;
    x->rpchop = 0.5;
    x->normalize = 0;
    x->skip = 0;
    x->bins = 0;
    x->logmag = 0;
    x->dbfloor = -120.0;
    x->state = new CccRpcState(static_cast<size_t>(highDim), static_cast<size_t>(maxLowDim),
                               static_cast<size_t>(maxFrame));

    // outlets are created right to left
    x->out_count = intout((t_object*)x);
    x->out_value = floatout((t_object*)x);

    attr_args_process(x, (short)argc, argv);
    return x;
}

void cccrpc_free(t_cccrpc* x) {
    object_free(x->bufref);
    delete x->state;
    x->state = nullptr;
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

t_max_err cccrpc_notify(t_cccrpc* x, t_symbol* s, t_symbol* msg, void* sender, void* data) {
    return x->bufref ? buffer_ref_notify(x->bufref, s, msg, sender, data) : MAX_ERR_NONE;
}

// Apply @skip / @bins / @logmag to state.raw, writing state.frame. Returns the length.
static size_t cccrpc_prepare(t_cccrpc* x, size_t n) {
    CccRpcState& st = *x->state;
    const size_t start = std::min(static_cast<size_t>(std::max<long>(0, x->skip)), n);
    size_t count = n - start;
    if (x->bins > 0) count = std::min(count, static_cast<size_t>(x->bins));
    count = std::min(count, st.maxFrame);

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
static double cccrpc_noise_reference(t_cccrpc* x, size_t n, size_t lowDim, size_t res, double hop) {
    CccRpcState& st = *x->state;
    if (st.refN == n && st.refLowDim == lowDim && st.refRes == res && st.refHop == hop) return st.noiseRef;

    cccrt::rpc::Pcg32 rng(1234);   // fixed seed: the same settings always give the same reference
    const int nFrames = 8;
    double sum = 0;
    for (int f = 0; f < nFrames; ++f) {
        for (size_t i = 0; i < n; ++i) st.calibFrame[i] = rng.uniform01<double>() * 2.0 - 1.0;
        sum += cccrt::rpc::calc(st.matrix.data(), lowDim, st.highDim, st.calibFrame.data(), n,
                                res, hop, st.calibProj.data(), st.calibCells.data());
    }
    st.noiseRef = std::max(1.0, sum / nFrames);
    st.refN = n; st.refLowDim = lowDim; st.refRes = res; st.refHop = hop;
    return st.noiseRef;
}

// Analyse state.frame[0..count) and send the result out.
static void cccrpc_analyse(t_cccrpc* x, size_t count) {
    CccRpcState& st = *x->state;
    outlet_int(x->out_count, static_cast<t_atom_long>(count));
    if (count < st.highDim) {   // too short to project even once
        outlet_float(x->out_value, 0.0);
        return;
    }
    const size_t lowDim = std::min(static_cast<size_t>(std::max<long>(1, x->lowdim)), st.maxLowDim);
    const size_t res = static_cast<size_t>(std::max<long>(1, x->res));
    const double hop = x->rpchop;

    double value = cccrt::rpc::calc(st.matrix.data(), lowDim, st.highDim, st.frame.data(), count,
                                    res, hop, st.proj.data(), st.cells.data());
    if (x->normalize == 1) {
        const double mx = cccrt::rpc::maxOccupiedCells(count, st.highDim, hop, res, lowDim);
        value = mx > 0 ? value / mx : 0;
    } else if (x->normalize == 2) {
        value /= cccrpc_noise_reference(x, count, lowDim, res, hop);
    }
    outlet_float(x->out_value, value);
}

void cccrpc_list(t_cccrpc* x, t_symbol* s, long argc, t_atom* argv) {
    CccRpcState& st = *x->state;
    const size_t n = std::min(static_cast<size_t>(std::max<long>(0, argc)), st.maxFrame);
    if (static_cast<size_t>(argc) > st.maxFrame) {
        object_warn((t_object*)x, "frame of %ld values truncated to maxFrame (%ld)", argc, (long)st.maxFrame);
    }
    for (size_t i = 0; i < n; ++i) st.raw[i] = atom_getfloat(argv + i);
    cccrpc_analyse(x, cccrpc_prepare(x, n));
}

void cccrpc_bang(t_cccrpc* x) {
    CccRpcState& st = *x->state;
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
    const size_t n = std::min(frames, st.maxFrame);
    for (size_t i = 0; i < n; ++i) st.raw[i] = samples[i * channels];   // channel 0
    buffer_unlocksamples(buffer);

    cccrpc_analyse(x, cccrpc_prepare(x, n));
}
