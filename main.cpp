#include <iostream>
#include <armadillo>
#include "gtest/gtest.h"
#include "shannonEntropy.hpp"
#include "ETC.hpp"
#include <thread>
#include "CCC.hpp"
#include <Eigen/Dense>
 
using Eigen::ArrayXi;

using namespace std;
using namespace arma;

// The fixture for testing class Foo.
//class CCCTest : public ::testing::Test {
//protected:
//    // You can remove any or all of the following functions if its body
//    // is empty.
//
//    CCCTest() {
//        // You can do set-up work for each test here.
//    }
//
//    ~CCCTest() override {
//        // You can do clean-up work that doesn't throw exceptions here.
//    }
//
//    // If the constructor and destructor are not enough for setting up
//    // and cleaning up each test, you can define the following methods:
//
//    void SetUp() override {
//        // Code here will be called immediately after the constructor (right
//        // before each test).
//    }
//
//    void TearDown() override {
//        // Code here will be called immediately after each test (right
//        // before the destructor).
//    }
//    // Objects declared here can be used by all tests in the test suite for Foo.
//};

ArrayXL ei(std::vector<long> v) {
    ArrayXL seq(v.size()); 
    for(size_t i=0; i < v.size(); i++) seq[i] = v[i];
    return seq;
}

TEST(CCCTest, ShannonEntropyTest) {

    EXPECT_DOUBLE_EQ(shannonEntropy::calc(ei({1,0,0,0,1,100,1,101,2,4,5,1})), 2.522055208874201);
    EXPECT_DOUBLE_EQ(shannonEntropy::calc(ei({1,0,0,0})), 0.811278124459133);
    EXPECT_DOUBLE_EQ(shannonEntropy::calc(ei({0})), 0);
    EXPECT_DOUBLE_EQ(shannonEntropy::calc(ei({1,0,3,1,3,1,4,0,2,9,4,1,2,4,2,2,34,4,100,300,-20,20,-111,3,0,2,0,0,1,4,0,19,34,235,7,27,25,43})), 3.745464778448826);
}

TEST(CCCTest, PairTest) {
    EXPECT_EQ(ETC::findHFPair(ei({1,0,1,0,3,3})).i128, ETC::makeETCPair(1,0).i128);
    EXPECT_EQ(ETC::findHFPair(ei({1,1,1,0,3,3,0,2,2,4,3,3})).i128, ETC::makeETCPair(3,3).i128);
    EXPECT_EQ(ETC::findHFPair(ei({1,1,1,1,0,3,3,0,2,2,4,3,3})).i128, ETC::makeETCPair(1,1).i128);
}

bool arma_eq(const ArrayXL &v1, const ArrayXL &v2) {
    bool res = 1;
    for(size_t i=0; i < v1.size(); i++) {
        if (v1[i] != v2[i])
        {
            res = 0;
            break;
        }
    }
    return res;
}

TEST(CCCTest, SubstituteTest) {
    EXPECT_EQ(arma_eq(get<0>(ETC::substitute(ei({1,1,3}), ETC::makeETCPair(1,1))), ei({4,3})), 1);
    EXPECT_EQ(arma_eq(get<0>(ETC::substitute(ei({-1,-1,-3}), ETC::makeETCPair(-1,-3))),ei({-1,0})), 1);
    auto [newSeqRepl, replaceCount, replaceSym]  = ETC::substitute(ei({1,1,3,3,5,8,6,4,3,6,3,6,2,3,1,3,5,3,6,4,2,3,4,5,1,3,5,4,2,3}), ETC::makeETCPair(1,3));
    EXPECT_EQ(arma_eq(newSeqRepl, 
        ei({1,9,3,5,8,6,4,3,6,3,6,2,3,9,5,3,6,4,2,3,4,5,9,5,4,2,3}))
        , 1);
    EXPECT_EQ(arma_eq(get<0>(ETC::substitute(ei({4,7,8,8,5,2,8,8}), ETC::makeETCPair(8,8))),ei({4,7,9,5,2,9})), 1);
}


TEST(CCCTest, ETCTest) {
    EXPECT_DOUBLE_EQ(ETC::calc(ei({1,0,1,0,3,3})),4/5.0);
    EXPECT_DOUBLE_EQ(ETC::calc(ei({1,4,5,3,3,3,5,5,7,5,7,7,7,7,8,2,4,9,2,4,6,7,8,5})),20/23.0);
    EXPECT_DOUBLE_EQ(ETC::calc(ei({1,4,5,3, 3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3})),8/22.0);
    EXPECT_DOUBLE_EQ(ETC::calc(ei({4294967297,21474836482,12884901891,17179869188,4294967297,8589934594})),1);
    
}
TEST(CCCTest, ETCJointTest) {
    EXPECT_DOUBLE_EQ(ETC::calcJoint(ei({1,1,1,2}), ei({2,3,3,4})),3/3.0);
    EXPECT_DOUBLE_EQ(ETC::calcJoint(ei({3,33,333,22, 30000,-20000,16,0}),ei({0,333,33,30000, 30000,30000,64,-1})),7/7.0);
    EXPECT_DOUBLE_EQ(ETC::calcJoint(ei({0,0,0,0,0}),ei({1,1,1,1,1})),0/4.0);
    EXPECT_DOUBLE_EQ(ETC::calcJoint(ei({0,1,2,2, 4,5,6,7, 2,2}),ei({0,1,2,2,4,5,6,7,2,2})),8/9.0);
    EXPECT_DOUBLE_EQ(ETC::calcJoint(ei({0}), ei({100})),0);
    
}

TEST(CCCTest, DCTest) {
    EXPECT_FLOAT_EQ(CCC::dynamicCC(ei({18,24,18,15,10,16,10,11,10,12,13,23,6,9,15,13,7,15,22,2,12,3,15,14,20,6,15,12,22,17,9,9,7,24,8,4,10,19,18,18,1,10,23,8,15,8,4,10,21,25,22,3,9,6,8,25,14,19,22,5,23,3,19,19,18,4,12,13,14,22,17,21,14,24,2,14,8,13,18,6,16,9,23,4,3,24,4,4,4,3,12,17,21,20,18,12,18,24,13,8,20,6,6,12,16,16,4,4,8,19,11,21,24,10,2,15,19,9,13,6,24,13,14,6,3,2,21,20,5,23,18,17,25,1,11,13,15,3,17,2,14,18,13,16,6,17,10,3,10,7,7,16,14,11,6,22,22,8,16,20,24,23,10,5,20,3,4,9,22,15,15,24,15,1,21,16,12,7,7,13,}), 7, 3, 1), 0.003267973856209149124963);
    EXPECT_FLOAT_EQ(CCC::dynamicCC(ei({7,10,9,9,8,10,9,3,11,2,8,3,3,4,3,6,11,5,7,7,8,5,10,11,9,4,9,2,6,7,8,5,9,10,1,2,10,5,5,11,9,11,3,2,5,6,7,11,6,5,11,3,8,11,8,4,6,1,10,5,4,11,6,4,1,7,2,1,7,4,}), 10, 17, 3), -0.045192307692307663591347);
    EXPECT_FLOAT_EQ(CCC::dynamicCC(ei({12,3,10,8,12,8,9,9,5,4,1,6,8,4,7,3,4,8,11,11,3,5,11,4,12,8,2,10,8,7,4,1,3,5,8,12,3,10,11,6,2,10,6,2,11,8,2,11,1,1,12,9,5,7,10,1,9,9,6,7,5,3,5,12,3,5,4,2,7,10,11,9,5,4,9,11,12,1,6,2,4,11,3,2,7,4,5,10,11,9,6,8,3,8,9,8,5,5,3,10,6,5,10,3,11,10,11,5,5,11,6,7,9,12,7,8,7,12,11,3,3,3,10,3,12,9,3,11,11,12,7,7,3,7,7,2,6,7,3,8,1,4,7,4,8,5,2,}), 21, 15, 2), -0.057908162);
    
    EXPECT_FLOAT_EQ(CCC::dynamicCC(ei({7,2,2,2,6,2,3,7,5,8,1,4,5,5,7,2,3,1,8,7,7,2,3,2,4,7,6,6,6,6,7,6,2,5,2,8,1,3,5,4,7,5,2,6,2,6,6,2,2,2,2,2,}), 22, 3, 3), -0.134259259259259272623055);
    EXPECT_FLOAT_EQ(CCC::dynamicCC(ei({2,4,6,6,2,4,6,6,1,2,5,1,1,3,6,6,1,4,6,6,6,6,5,1,4,2,5,6,2,6,3,3,2,6,5,3,5,6,2,2,4,5,6,3,1,1,3,5,2,5,4,1,5,6,3,6,1,5,6,5,2,2,4,3,3,6,4,6,3,3,5,5,3,1,1,5,4,1,4,5,1,2,5,1,6,6,1,3,2,4,1,1,6,6,3,4,1,2,1,}), 12, 6, 3), -0.099782132);
    
}

TEST(CCCTest, DCJointTest) {
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({1,2,3,4,1,2,3,4,1,2,3,4}), ei({1,5,3,7,1,1,3,7,1,2,4,7}), 3, 3, 1), 0);
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({1,2,3,4,1,1,3,7,1,2,3,7}), ei({1,5,3,7,1,1,3,7,1,2,4,7}), 3, 3, 1), 0);
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({1,2,3,4,1,1,3,7,1,2,3,7}), ei({1,5,3,7,1,1,3,7,1,2,4,7}), 4, 4, 1), 0);
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({4,4,1,1,1,2,1,2,2,3,4,2,2,4,2,1,3,4,2,3,3,2,2,4,1,2,2,3,4,2,2,3,2,4,2,2,2,1,1,1,4,1,4,3,4,4,1,2,4,3,}), ei({4,2,1,3,1,3,3,3,3,2,1,2,1,3,4,1,1,4,1,1,4,1,4,1,2,3,3,3,3,2,1,4,1,3,4,1,4,1,4,4,1,3,4,2,4,3,4,4,1,2,}), 9, 12, 3), -0.050000000000000023592239);
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({3,3,2,2,4,2,4,1,4,1,4,4,3,1,1,3,4,1,4,4,2,3,2,1,4,3,2,3,2,1,1,2,2,2,2,1,2,3,1,3,2,2,1,1,3,2,3,4,3,3,}), ei({4,4,3,4,1,2,4,4,2,3,2,2,1,2,1,2,2,4,3,4,2,2,4,4,3,4,3,3,4,4,4,3,2,2,4,1,1,2,3,1,3,4,1,2,2,4,3,2,1,3,}), 11, 19, 3), -0.048987411056376557738634);
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({1,3,4,1,4,4,1,1,2,4,2,2,3,1,1,4,3,2,4,1,2,3,2,2,3,4,4,1,2,2,4,2,1,3,2,4,4,4,3,3,4,1,2,1,4,3,2,1,3,2,}), ei({1,1,1,2,1,2,4,2,1,3,4,2,1,1,1,3,2,3,3,3,3,4,3,4,2,1,1,3,2,3,3,3,3,3,2,2,4,2,1,1,3,4,4,1,4,3,2,3,2,3,}), 10, 9, 1), -0.035842293906810054893164);
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({4,2,4,1,1,4,4,4,4,1,3,4,3,4,1,4,4,3,2,2,3,1,2,4,3,2,2,2,4,4,3,4,2,2,1,1,4,3,3,4,3,3,1,3,1,1,3,3,1,2,}), ei({3,3,1,1,1,2,2,1,3,1,1,4,4,1,2,4,1,4,4,3,1,2,4,3,2,3,1,3,4,4,2,3,2,4,1,2,4,1,3,4,3,1,2,2,2,3,1,4,3,1,}), 8, 21, 3), -0.030612244897959169087631);
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({3,2,2,3,1,2,4,4,4,4,1,3,3,4,3,2,3,4,3,1,4,1,3,1,2,4,1,4,3,2,1,2,4,1,3,3,3,4,2,2,4,4,2,1,3,1,3,2,2,3,}), ei({3,4,1,3,1,4,3,4,1,4,3,3,2,1,4,1,2,3,4,3,4,1,2,3,1,3,2,4,2,2,4,3,2,1,3,3,1,1,4,2,1,1,1,4,2,3,2,1,4,2,}), 5, 14, 2), -0.031250000000000013877788);
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({2,3,4,3,2,1,2,1,4,3,1,3,1,4,2,3,3,3,2,1,2,1,3,1,4,3,4,3,1,2,4,4,2,4,1,3,2,3,2,1,2,3,1,2,3,4,2,2,4,3,}), ei({1,4,2,2,1,2,2,2,2,1,4,4,2,2,3,2,2,4,1,4,1,2,3,3,1,4,3,2,3,4,4,2,2,2,1,4,1,2,2,4,1,4,2,2,2,3,1,1,4,3,}), 6, 22, 3), -0.010582010582010595300950);
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({2,1,5,7,1,1,5,5,4,1,5,7,7,1,3,2,4,5,3,6,6,6,3,1,5,7,2,1,7,4,6,6,2,2,1,5,6,4,1,7,4,3,2,5,5,6,2,3,2,6,}), ei({3,2,1,6,4,4,7,6,5,2,4,7,5,6,7,7,2,1,4,1,4,7,2,3,4,6,6,4,3,7,6,6,7,4,6,4,6,3,6,7,4,4,2,3,1,3,3,3,6,4,}), 8, 18, 1), -0.003333333333333336149368);
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint(ei({2,2,2,2,3,1,3,3,1,1,2,3,2,1,2,2,2,1,3,2,2,1,2,1,1,3,2,1,2,3,1,3,1,1,1,3,1,3,1,1,3,1,2,3,3,3,3,2,1,2,}), ei({1,3,2,3,3,2,1,3,3,3,3,3,2,3,2,3,1,3,1,2,2,1,3,2,3,2,3,2,3,3,2,3,2,3,1,3,1,1,3,3,3,1,3,2,2,3,3,1,3,1,}), 6, 8, 3), -0.078754578754578710708678);
    
}


TEST(CCCTest, CCCausalityTest) {
// //    EXPECT_FLOAT_EQ(CCC::CCCausality({1,2,3,4,1,1,3,7,1,2,3,7,7,}, {1,5,3,7,1,1,3,7,1,2,4,7,3,}, 3, 3, 1), 0);
    EXPECT_FLOAT_EQ(get<0>(CCC::CCCausality(ei({1,5,3,7,1,1,3,7,1,2,4,7,3,}),ei({1,2,3,4,1,1,3,7,1,2,3,7,7,}), 3, 3, 1)), -0.057142857142857127195068);
    EXPECT_FLOAT_EQ(get<0>(CCC::CCCausality(ei({2,4,3,3,3,3,1,4,4,1,4,4,4,1,4,3,1,2,4,4,1,2,1,4,2,4,2,4,1,4,1,3,3,2,1,4,1,2,1,1,2,4,1,3,3,1,2,4,1,3,}), ei({3,3,2,3,2,1,2,3,3,4,4,1,3,3,2,2,2,3,4,1,4,3,4,1,2,1,3,4,2,1,2,2,3,3,2,4,4,3,1,3,2,4,2,3,3,1,3,1,1,1,}), 3, 3, 1)), -0.009090909090909073120290);
    EXPECT_FLOAT_EQ(get<0>(CCC::CCCausality(ei({1,2,3,2,4,2,1,3,1,4,4,2,1,1,2,4,3,2,2,2,4,4,3,3,3,1,3,1,2,2,1,1,4,3,2,2,1,1,2,2,4,2,2,1,4,3,4,3,4,4,}), ei({2,3,3,2,4,1,1,1,1,3,3,2,2,1,4,3,2,1,4,4,4,2,4,4,1,1,4,2,2,4,4,4,3,4,2,1,2,2,3,1,1,1,4,4,3,4,4,4,1,4,}), 3, 3, 1)), 0.009090909090909099141142);
    EXPECT_FLOAT_EQ(get<0>(CCC::CCCausality(ei({1,3,4,2,4,1,2,2,1,3,2,4,2,3,1,3,3,4,3,4,4,1,1,4,1,2,2,3,1,4,4,3,2,1,3,1,1,1,2,3,4,2,1,2,3,4,2,2,4,3,}), ei({3,2,2,4,2,2,2,3,3,3,2,2,2,1,2,4,4,3,4,4,4,4,4,1,2,2,2,1,4,3,1,1,3,2,1,2,1,4,4,1,3,3,2,4,2,1,1,4,3,4,}), 9, 10, 1)), -0.014336917562723941466096);
    EXPECT_FLOAT_EQ(get<0>(CCC::CCCausality(ei({3,4,3,4,4,4,1,1,1,1,1,3,3,4,3,2,2,2,3,1,3,4,3,2,3,2,1,4,3,2,4,1,2,2,1,1,3,1,3,2,2,4,1,1,4,3,3,3,1,4,}), ei({2,3,2,1,3,3,1,1,4,1,2,4,4,3,4,3,3,3,4,3,2,3,2,1,2,3,2,2,4,3,2,1,2,2,3,3,4,4,1,4,1,4,3,1,1,2,4,1,2,1,}), 6, 6, 3)), -0.020979020979020986809038);
    EXPECT_FLOAT_EQ(get<0>(CCC::CCCausality(ei({2,2,2,1,1,2,1,2,1,1,1,2,2,2,2,2,2,2,1,1,1,2,2,1,2,1,1,2,2,2,2,2,2,1,1,2,1,1,1,2,2,2,2,2,2,1,2,1,1,1,}), ei({2,2,2,1,2,2,1,1,2,2,1,1,1,2,1,2,2,1,2,2,2,2,1,1,1,1,1,1,2,1,2,2,1,1,1,2,1,2,1,1,2,2,2,1,1,2,2,1,1,1,}), 8, 6, 1)), 0.027777778);

    EXPECT_FLOAT_EQ(get<0>(CCC::CCCausality(ei({10,7,7,10,5,1,14,4,17,9,7,15,15,4,3,4,13,17,11,4,15,19,5,18,5,17,18,11,10,18,2,11,13,15,5,1,3,16,20,8,20,7,3,8,12,13,11,7,4,4,}), ei({3,15,10,3,15,4,17,11,11,5,16,11,12,9,15,2,12,18,9,14,7,13,6,16,16,20,9,10,12,17,1,4,19,19,10,17,19,14,8,9,9,14,10,20,8,7,18,11,15,3,}), 8, 6, 2)), 0);
    
    
}
TEST(CCCTest, CCCausalityMassTest) {
    ArrayXL s1(100),s2(100);
    for(int i=0; i < 10; i++) {
        for(int j=0; j < s1.size(); j++) {
            s1[j] = rand()  % 16;
            s2[j] = rand()  % 16;
        }
        double res = get<0>(CCC::CCCausality(s1, s2, 30, 50, 1));
        cout << res << ", ";
    }
}

int main(int argc, char **argv) {
    cout << "CCC library tests\n";
    // VectorXi test {0};
    ::testing::InitGoogleTest(&argc, argv);
    int res =  RUN_ALL_TESTS();
    unsigned int n = std::thread::hardware_concurrency();
    std::cout << n << " concurrent threads are supported.\n";
    //perf testing
//    clock_t t = clock();
//    for(int i=0; i < 86; i++) {
//        ivec x(86);
//        for(int j=0; j<x.size(); j++) x[j] = rand() % 16;
//        cout << ETC::calc(x) << endl;
//    }
//    const double work_time = (clock() - t) / double(CLOCKS_PER_SEC) * 1000;
//    cout << work_time << " ms" << endl;
    
//    ivec test = {0,1,2,3,4,5,6,7,8,9};
//    int offset=9, dx=2, past=3;
//    cout << test.subvec(offset-dx-past+1,offset) << endl;
//    cout << test.subvec(offset-dx-past+1,offset-dx) << endl;
//    ivec test = {1,2,3,4,5,6,7,8,9};
//    ivec test2 = {10,11,12,13,14,15};
//    ivec comb = test.subvec(0,6);
//    comb(span(4,5)) = test2(span(4,5));
//    cout << test.t() << endl;
//    cout << test2.t() << endl;
//    cout << comb.t() << endl;
    return res;
}


/*
 realtime stuff -
 spectral RMS rate, 512 hop size = 86Hz
 Looking at a window of 1 sec = 86 values, 86 win size, sym size =16 -> 49ms for worst case random data -->>> realtime is fine
 for CCC measure, can use large scale multithreading
 
 */
