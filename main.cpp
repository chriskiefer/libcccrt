#include <iostream>
#include <armadillo>
#include "gtest/gtest.h"
#include "shannonEntropy.hpp"
#include "ETC.hpp"
#include <thread>
#include "CCC.hpp"

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

TEST(CCCTest, ShannonEntropyTest) {
    EXPECT_DOUBLE_EQ(shannonEntropy::calc({1,0,0,0,1,100,1,101,2,4,5,1}), 2.522055208874201);
    EXPECT_DOUBLE_EQ(shannonEntropy::calc({1,0,0,0}), 0.811278124459133);
    EXPECT_DOUBLE_EQ(shannonEntropy::calc({0}), 0);
    EXPECT_DOUBLE_EQ(shannonEntropy::calc({1,0,3,1,3,1,4,0,2,9,4,1,2,4,2,2,34,4,100,300,-20,20,-111,3,0,2,0,0,1,4,0,19,34,235,7,27,25,43}), 3.745464778448826);
}

TEST(CCCTest, PairTest) {
    EXPECT_EQ(ETC::findHFPair({1,0,1,0,3,3}).i128, ETC::makeETCPair(1,0).i128);
    EXPECT_EQ(ETC::findHFPair({1,1,1,0,3,3,0,2,2,4,3,3}).i128, ETC::makeETCPair(3,3).i128);
    EXPECT_EQ(ETC::findHFPair({1,1,1,1,0,3,3,0,2,2,4,3,3}).i128, ETC::makeETCPair(1,1).i128);
}

bool arma_eq(const ivec &v1, const ivec &v2) {
    bool res = 1;
    for(uword i=0; i < v1.size(); i++) {
        if (v1[i] != v2[i])
        {
            res = 0;
            break;
        }
    }
    return res;
}

TEST(CCCTest, SubstituteTest) {
    EXPECT_EQ(arma_eq(get<0>(ETC::substitute({1,1,3}, ETC::makeETCPair(1,1))), {4,3}), 1);
    EXPECT_EQ(arma_eq(get<0>(ETC::substitute({-1,-1,-3}, ETC::makeETCPair(-1,-3))), {-1,0}), 1);
    EXPECT_EQ(arma_eq(get<0>(ETC::substitute({1,1,3,3,5,8,6,4,3,6,3,6,2,3,1,3,5,3,6,4,2,3,4,5,1,3,5,4,2,3}, ETC::makeETCPair(1,3))), {1,9,3,5,8,6,4,3,6,3,6,2,3,9,5,3,6,4,2,3,4,5,9,5,4,2,3}), 1);
}
TEST(CCCTest, ETCTest) {
    EXPECT_DOUBLE_EQ(ETC::calc({1,0,1,0,3,3}),4/5.0);
    EXPECT_DOUBLE_EQ(ETC::calc({1,4,5,3,3,3,5,5,7,5,7,7,7,7,8,2,4,9,2,4,6,7,8,5}),20/23.0);
    EXPECT_DOUBLE_EQ(ETC::calc({1,4,5,3, 3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3}),8/22.0);
}
TEST(CCCTest, ETCJointTest) {
    EXPECT_DOUBLE_EQ(ETC::calcJoint({1,1,1,2}, {2,3,3,4}),3/3.0);
    EXPECT_DOUBLE_EQ(ETC::calcJoint({3,33,333,22, 30000,-20000,16,0},{0,333,33,30000, 30000,30000,64,-1}),7/7.0);
    EXPECT_DOUBLE_EQ(ETC::calcJoint({0,0,0,0,0},{1,1,1,1,1}),0/4.0);
    EXPECT_DOUBLE_EQ(ETC::calcJoint({0,1,2,2, 4,5,6,7, 2,2},{0,1,2,2,4,5,6,7,2,2}),8/9.0);
    EXPECT_DOUBLE_EQ(ETC::calcJoint({0}, {100}),0);
}

TEST(CCCTest, DCTest) {
    EXPECT_FLOAT_EQ(CCC::dynamicCC({4,1,3,2,3,3,2,3,1,3,2,3,2,3,2,3,4,1,3,1}, 5, 12, 1), 0.073863636363636);
    EXPECT_FLOAT_EQ(CCC::dynamicCC({4,1,3,2,3,3,2,3,1,3,2,3,2,3,2,3,4,1,3,1}, 4, 13, 1), 0.020833333333333);
    EXPECT_FLOAT_EQ(CCC::dynamicCC({4,1,3,2,3,3,2,3,1,3,2,3,2,3,2,3,4,1,3,1}, 3, 13, 1), 0.016666666666667);
    EXPECT_FLOAT_EQ(CCC::dynamicCC({4,9,6,4,10,9,6,7,6,3}, 4, 4, 1), -0.285714285714286);
    EXPECT_FLOAT_EQ(CCC::dynamicCC({10,5,2,10,10,5,2,3,5,6,}, 3, 3, 1), -0.133333333333333);
    EXPECT_FLOAT_EQ(CCC::dynamicCC({4,2,5,5,2,6,3,4,6,3}, 2, 2, 2), 0.5);
    EXPECT_FLOAT_EQ(CCC::dynamicCC({3,7,3,9,10,8,4,6,2,10,9,9,3,6,1,5,4,2,2,5,1,6,5,7,7,7,1,1,4,6}, 2, 2, 2), 0.166666666666667);
    EXPECT_FLOAT_EQ(CCC::dynamicCC({6,3,7,1,7,7,8,9,10,8,6,10,6,1,2,9,5,9,3,6,7,1,7,4,1,5,2,2,3,2,}, 2, 2, 2), 8.333333e-02);
    EXPECT_FLOAT_EQ(CCC::dynamicCC({9,2,13,3,4,11,7,10,6,4,1,9,6,6,8,1,5,11,10,2,2,2,1,6,9,10,7,2,9,2,2,2,2,3,3,5,5,3,}, 5, 10, 3), -0.030612244897959186434866);

    EXPECT_FLOAT_EQ(CCC::dynamicCC({4,2,16,13,10,6,3,11,17,3,5,7,2,12,7,17,7,11,3,7,3,13,15,6,12,6,10,15,11,6,6,8,8,7,10,13,8,8,3,1,5,6,12,17,16,8,5,13,13,13,13,2,12,8,4,2,15,3,3,12,16,9,12,3,17,10,12,1,14,13,3,9,6,10,7,8,4,5,1,16,12,16,3,16,14,10,8,}, 8, 18, 1), -0.011529411764705883511328);
    
    EXPECT_FLOAT_EQ(CCC::dynamicCC({17,10,22,22,20,9,11,6,18,20,21,13,14,4,20,10,5,20,17,20,7,15,15,3,9,7,16,7,20,19,9,11,16,19,14,13,8,11,16,20,16,1,15,10,10,3,18,8,6,8,9,13,13,9,9,12,15,21,16,9,19,3,2,2,4,8,7,1,12,3,4,14,19,22,13,22,13,12,8,10,11,2,20,2,10,19,9,14,19,20,21,5,6,20,14,12,14,19,12,5,10,10,22,14,16,16,8,12,13,4,13,16,10,19,17,8,10,9,18,17,10,16,21,18,16,3,9,13,11,2,6,19,1,20,2,15,12,5,13,3,15,14,2,2,4,1,10,}, 19, 15, 2), -0.013914656771799618212304);
    
    EXPECT_NEAR(CCC::dynamicCC({2,1,7,4,7,6,1,5,5,5,4,2,6,2,3,7,6,3,3,5,7,7,5,3,6,4,7,1,4,6,2,3,2,3,3,4,1,4,2,2,2,2,7,7,6,6,2,3,2,1,3,5,5,4,4,3,4,6,6,5,6,5,1,4,3,1,2,2,5,4,5,2,1,4,2,7,7,3,1,5,6,7,1,4,5,5,6,5,6,3,5,1,1,7,2,5,7,2,2,3,7,3,2,2,3,3,1,4,1,5,1,5,6,2,4,4,1,4,5,2,6,7,6,4,2,6,2,7,5,5,2,1,2,7,7,5,6,2,}, 14, 19, 2), -0.054553952991453005805234, 0.001);


}

TEST(CCCTest, DCJointTest) {
    EXPECT_FLOAT_EQ(CCC::dynamicCCJoint({1,2,3,4,1,2,3,4,1,2,3,4}, {1,5,3,7,1,1,3,7,1,2,4,7}, 3, 3, 1), -0.04);
}

int main(int argc, char **argv) {
    cout << "CCC library tests\n";
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
