//
//  CCCpy.cpp
//  CCC
//
//  Created by Chris Kiefer on 06/10/2019.
//


#include <stdio.h>
#include "ETC.hpp"
#include "CCC.hpp"
#include "LZ.hpp"
#include "RPC.hpp"
#include "shannonEntropy.hpp"
#include "fractal.hpp"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

// namespace np = boost::python::numpy;
// using namespace boost::python;

namespace py = pybind11;
using namespace std;


char const* greet()
{
   return "Cheese!  It's something that grows on cheese trees?";
}

// double ETC(np::ndarray const & ar) {
//   int64_t* input_ptr = reinterpret_cast<int64_t*>(ar.get_data());
//   ivec inputVec(input_ptr, ar.shape(0), 1, 1);
//   return ETC::calc(inputVec);
// }

// double dynamicCC(np::ndarray const & ar, size_t dx, size_t xpast, size_t step) {
//   int64_t* input_ptr = reinterpret_cast<int64_t*>(ar.get_data());
//   ivec inputVec(input_ptr, ar.shape(0), 1, 1);
//   return CCC::dynamicCC(inputVec, dx, xpast, step);
// }

// double dynamicCCJoint(np::ndarray const & ar, np::ndarray const & ar2, size_t dx, size_t past, size_t step) {
//   int64_t* input_ptr = reinterpret_cast<int64_t*>(ar.get_data());
//   ivec inputVec(input_ptr, ar.shape(0), 1, 1);
//   int64_t* input_ptr2 = reinterpret_cast<int64_t*>(ar2.get_data());
//   ivec inputVec2(input_ptr2, ar.shape(0), 1, 1);
//   return CCC::dynamicCCJoint(inputVec, inputVec2, dx, past, step);
// }

// boost::python::tuple CCCausality(np::ndarray const & ar, np::ndarray const & ar2, size_t dx, size_t past, size_t step) {
//   int64_t* input_ptr = reinterpret_cast<int64_t*>(ar.get_data());
//   ivec inputVec(input_ptr, ar.shape(0), 1, 1);

//   int64_t* input_ptr2 = reinterpret_cast<int64_t*>(ar2.get_data());
//   ivec inputVec2(input_ptr2, ar2.shape(0), 1, 1);

//   auto CCC = CCC::CCCausality(inputVec, inputVec2, dx, past, step);
//   return boost::python::make_tuple(get<0>(CCC), get<1>(CCC));
// }

double ETC(Eigen::Ref<ArrayXL> v) {
  return ETC::calc(v);
}

double ETCJoint(const Eigen::Ref<ArrayXL> v0, const Eigen::Ref<ArrayXL> v1) {
  return ETC::calcJoint(v0, v1);
}

double dynamicCC(const Eigen::Ref<ArrayXL> seq, const size_t dx, const size_t xpast, const size_t step) {
  return CCC::dynamicCC(seq, dx, xpast, step);
}

double dynamicCCJoint(const Eigen::Ref<ArrayXL> X, const Eigen::Ref<ArrayXL> Y, const size_t dx, const size_t past, const size_t step) {
  return CCC::dynamicCCJoint(X, Y, dx, past, step, CCC::SINGLETHREAD);
}

std::tuple<double, unsigned int> CCCausality(const Eigen::Ref<ArrayXL>  effectSeq, const Eigen::Ref<ArrayXL> causeSeq, const size_t dx, const size_t past, const size_t step) {
  return CCC::CCCausality(effectSeq, causeSeq, dx, past, step);
}

size_t LZ(const Eigen::Ref<ArrayXL> v) {
  return LZ::calc(v);
}

double NLZ(const Eigen::Ref<ArrayXL> v) {
  return LZ::calcNorm(v);
}

double shannonEntropy(const Eigen::Ref<ArrayXL> v) {
  return shannonEntropy::calc(v);
  
}

double randomProjectionComplexity(const Eigen::Ref<Eigen::MatrixXd> &projectionMatrix, const Eigen::Ref<Eigen::VectorXd> &data, const size_t resolution, const double hop=0.5) {
  return RPC::calc(projectionMatrix, data, resolution, hop);
  // return 0;
}

double sevcik(const Eigen::Ref<Eigen::VectorXd> v) {
  return fractal::sevcik::calc(v);
  
}

PYBIND11_MODULE(cccrt, m) {
    m.doc() = "Libcccrt"; // optional module docstring
    m.def("greet", &greet, "Hello");
    m.def("ETC", &ETC, "Calculates Effort To Compress on an array of symbols");
    m.def("ETCJoint", &ETCJoint, "Calculates joint Effort To Compress on two arrays of symbols");
    m.def("dynamicCC", &dynamicCC, "Calculates dynamical complexity on an array of symbols");
    m.def("dynamicCCJoint", &dynamicCCJoint, "Calculates joint dynamical complexity on an two arrays of symbols");
    m.def("CCCausality", &CCCausality, "Calculates joint compression-complexity causality on an two arrays of symbols");
    m.def("LZ", &LZ, "Calculates lempel-ziv complexity on an array of symbols");
    m.def("NLZ", &NLZ, "Calculates normalised lempel-ziv complexity on an array of symbols");
    m.def("shannonEntropy", &shannonEntropy, "Calculates shannon entropy on an array of symbols");
    m.def("createRPMatrix", &RPC::createProjectionMatrix, py::return_value_policy::reference_internal);
    m.def("RPC", &randomProjectionComplexity, "Random projection complexity");
    m.def("sevcik", &sevcik, "Sevcik fractal complexity");
}
