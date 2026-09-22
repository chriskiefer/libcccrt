// Standalone tests for the Eigen-free core.
// Builds with nothing but a C++17 compiler; no gtest, no Eigen, no heap use in the algorithms.
// The expected values are the same ones used in the main gtest suite (../../main.cpp).
#include "../cccrt.hpp"
#include <cstdio>
#include <cstdint>
#include <cmath>

static int failures = 0;

#define CHECK_NEAR(expr, expected, tol)                                                       \
    do {                                                                                      \
        const double got_ = (expr);                                                           \
        if (!(std::fabs(got_ - (expected)) <= (tol))) {                                       \
            std::printf("FAIL %s:%d: %s = %.15g, expected %.15g\n", __FILE__, __LINE__, #expr, \
                        got_, (double)(expected));                                            \
            ++failures;                                                                       \
        }                                                                                     \
    } while (0)

#define CHECK_EQ(expr, expected) CHECK_NEAR(expr, expected, 0)

static void testShannon() {
    int64_t a[] = {1, 0, 0, 0, 1, 100, 1, 101, 2, 4, 5, 1};
    int64_t scratch[64];
    CHECK_NEAR(cccrt::shannonEntropy(a, 12, scratch), 2.522055208874201, 1e-12);
    int64_t b[] = {1, 0, 0, 0};
    CHECK_NEAR(cccrt::shannonEntropy(b, 4, scratch), 0.811278124459133, 1e-12);
    int64_t c[] = {0};
    CHECK_EQ(cccrt::shannonEntropy(c, 1, scratch), 0);
    int64_t d[] = {1, 0, 3, 1, 3, 1, 4, 0, 2, 9, 4, 1, 2, 4, 2, 2, 34, 4, 100, 300, -20, 20, -111, 3, 0, 2, 0, 0, 1, 4, 0, 19, 34, 235, 7, 27, 25, 43};
    CHECK_NEAR(cccrt::shannonEntropy(d, 38, scratch), 3.745464778448826, 1e-12);
    // float scalar, small symbol type, in place
    uint8_t e[] = {1, 0, 0, 0};
    CHECK_NEAR(cccrt::shannonEntropyInPlace<float>(e, 4), 0.811278124459133, 1e-6);
}

static void testLZ() {
    uint32_t starts[64];
    int64_t a[] = {0, 1};
    CHECK_EQ(cccrt::lempelZiv(a, 2, starts), 2);
    int64_t b[] = {0, 1, 3, 3};
    CHECK_EQ(cccrt::lempelZiv(b, 4, starts), 3);
    int64_t c[] = {1, 2, 3, 3, 3, 2, 3, 1};
    CHECK_EQ(cccrt::lempelZiv(c, 8, starts), 5);
    int64_t d[] = {100, 9, 8, 101, 102, 2, 3, 4, 9, 100, 103, 104, 101, 105, 8, 9, 106, 2, 3, 4, 9, 8, 6, 105, 3, 100};
    CHECK_EQ(cccrt::lempelZiv(d, 26, starts), 19);
    uint8_t e[] = {1, 2, 3, 3, 3, 2, 3, 1};
    CHECK_EQ(cccrt::lempelZiv(e, 8, starts), 5);
    CHECK_NEAR(cccrt::lempelZivNorm<float>(e, 8, starts), 5.0 / (8.0 / std::log(8.0)), 1e-6);
}

static void testSevcik() {
    double a[] = {1.0, 3.0, 4.2, 4.2, 1.3, 1.2, 9.7, 4.2, 2.71, 2, 5, 10.3};
    CHECK_NEAR(cccrt::sevcik(a, 12), 1.417138504999742, 1e-12);
    double b[] = {0.59268124, 0.10992054, 0.49541226, 0.92321486, 0.15333986,
                  0.10439158, 0.42802243, 0.22755213, 0.71985698, 0.84432153};
    CHECK_NEAR(cccrt::sevcik(b, 10), 1.4940149591892897, 1e-5);
    float bf[10];
    for (int i = 0; i < 10; ++i) bf[i] = (float)b[i];
    CHECK_NEAR(cccrt::sevcik(bf, 10), 1.4940149591892897, 1e-5);
    float sine[100];
    for (int i = 0; i < 100; ++i) sine[i] = std::sin(i * 0.05f);
    CHECK_NEAR(cccrt::sevcik(sine, 100), 1.1202632289010772, 1e-5);
}

static void testRingBuffer() {
    float storage[8];
    cccrt::RingBuffer<float> ring(storage, 8);
    float out[8];
    // before anything is pushed, the buffer reads as zeros
    ring.copyLatest(out, 3);
    CHECK_EQ(out[0] + out[1] + out[2], 0);
    for (int i = 1; i <= 5; ++i) ring.push((float)i);
    ring.copyLatest(out, 3);            // most recent 3: 3 4 5
    CHECK_EQ(out[0], 3); CHECK_EQ(out[1], 4); CHECK_EQ(out[2], 5);
    ring.copyLatest(out, 2, 1);         // skipping the last one: 3 4
    CHECK_EQ(out[0], 3); CHECK_EQ(out[1], 4);
    for (int i = 6; i <= 13; ++i) ring.push((float)i);   // wraps around
    ring.copyLatest(out, 7);            // 7..13
    for (int i = 0; i < 7; ++i) CHECK_EQ(out[i], 7 + i);
    ring.copyLatest(out, 4, 3);         // 7..10
    for (int i = 0; i < 4; ++i) CHECK_EQ(out[i], 7 + i);
}

static void testRPC() {
    // flat index
    double i1[] = {3};
    CHECK_EQ(cccrt::rpc::flatIndex(i1, 1, 10), 3);
    double i2[] = {9, 9};
    CHECK_EQ(cccrt::rpc::flatIndex(i2, 2, 10), 99);
    double i3[] = {2, 2, 1};
    CHECK_EQ(cccrt::rpc::flatIndex(i3, 3, 4), 41);
    int i3b[] = {1, 1, 1};
    CHECK_EQ(cccrt::rpc::flatIndex(i3b, 3, 4), 21);

    // matrix statistics: mean ~0, sd ~ 1/sqrt(nDim)
    static float m[4 * 2000];
    cccrt::rpc::makeProjectionMatrix(m, 4, 2000);
    double sum = 0, sumsq = 0;
    for (int i = 0; i < 8000; ++i) { sum += m[i]; sumsq += m[i] * m[i]; }
    CHECK_NEAR(sum / 8000, 0.0, 0.02);
    CHECK_NEAR(std::sqrt(sumsq / 8000), 0.5, 0.02);
    // deterministic for a given seed
    static float m2[4 * 2000];
    cccrt::rpc::makeProjectionMatrix(m2, 4, 2000);
    CHECK_EQ(m[123], m2[123]);

    // normalisation bound: min(hops, res^dims)
    CHECK_EQ(cccrt::rpc::maxOccupiedCells(100, 4, 0.5, 5, 2), 25);     // 49 hops, 25 cells
    CHECK_EQ(cccrt::rpc::maxOccupiedCells(100, 4, 0.5, 100, 1), 49);   // 49 hops, 100 cells
    CHECK_EQ(cccrt::rpc::maxOccupiedCells(3, 4, 0.5, 5, 2), 0);        // shorter than the window

    // same bounds as RPCProjectionTest in main.cpp
    cccrt::rpc::Fixed<float, 3, 16, 128> rpc;
    if (!rpc.init(2, 4)) { std::puts("FAIL: rpc.init"); ++failures; }
    float sparse[] = {0, 0, 1, 0, 2, 0, 0, 0, 1, 0, 0, 0};
    const float lowComplexity = rpc.calc(sparse, 12, 10);
    if (!(lowComplexity < 5)) { std::printf("FAIL: rpc on sparse data = %g, expected < 5\n", lowComplexity); ++failures; }

    float noise[100];
    cccrt::rpc::Pcg32 rng(99);
    for (auto& v : noise) v = rng.uniform01<float>() * 2 - 1;
    const float highComplexity = rpc.calc(noise, 100, 4);
    if (!(highComplexity > 5)) { std::printf("FAIL: rpc on noise = %g, expected > 5\n", highComplexity); ++failures; }

    cccrt::rpc::Fixed<float, 3, 16, 128> rpc3;
    rpc3.init(3, 10);
    const float highComplexity3 = rpc3.calc(noise, 100, 3, 0.001f);
    if (!(highComplexity3 > 5)) { std::printf("FAIL: rpc3 on noise = %g, expected > 5\n", highComplexity3); ++failures; }

    // capacity guards
    static float longNoise[1000];
    for (auto& v : longNoise) v = rng.uniform01<float>() * 2 - 1;
    CHECK_EQ(rpc3.calc(longNoise, 1000, 3, 0.001f), 0);  // 991 hops > MaxHops -> 0
    CHECK_EQ(rpc.calc(sparse, 3, 10), 0);                 // shorter than window -> 0
    cccrt::rpc::Fixed<float, 2, 8, 8> tooSmall;
    if (tooSmall.init(3, 4)) { std::puts("FAIL: init should reject nDim > MaxDims"); ++failures; }
}

int main() {
    testShannon();
    testLZ();
    testSevcik();
    testRingBuffer();
    testRPC();
    if (failures == 0) std::puts("core tests: all passed");
    else std::printf("core tests: %d failure(s)\n", failures);
    return failures != 0;
}
