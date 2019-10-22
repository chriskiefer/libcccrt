#include <iostream>
#include <armadillo>
#include "gtest/gtest.h"
#include "shannonEntropy.hpp"
#include "ETC.hpp"

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
    EXPECT_EQ(ETC::calc({1,0,1,0,3,3}),4);
    EXPECT_EQ(ETC::calc({1,4,5,3,3,3,5,5,7,5,7,7,7,7,8,2,4,9,2,4,6,7,8,5}),20);
    EXPECT_EQ(ETC::calc({1,4,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3}),8);
}
TEST(CCCTest, ETCJointTest) {
    EXPECT_EQ(ETC::calcJoint({1,1,1,2}, {2,3,3,4}),3);
    EXPECT_EQ(ETC::calcJoint({3,33,333,22,30000,-20000,16,0},{0,333,33,30000,30000,30000,64,-1}),7);
    EXPECT_EQ(ETC::calcJoint({0,0,0,0,0},{1,1,1,1,1}),0);
    EXPECT_EQ(ETC::calcJoint({0,1,2,2,4,5,6,7,2,2},{0,1,2,2,4,5,6,7,2,2}),8);
    EXPECT_EQ(ETC::calcJoint({0}, {100}),0);


}


int main(int argc, char **argv) {
    cout << "CCC library tests\n";
    ::testing::InitGoogleTest(&argc, argv);
    int res =  RUN_ALL_TESTS();
    
    //perf testing
//    clock_t t = clock();
//    for(int i=0; i < 50; i++) {
//        ivec x(500);
//        for(int j=0; j<x.size(); j++) x[j] = rand() % 8;
//        cout << ETC::calc(x) << endl;
//    }
//    const double work_time = (clock() - t) / double(CLOCKS_PER_SEC) * 1000;
//    cout << work_time << " ms" << endl;
    return res;
}
