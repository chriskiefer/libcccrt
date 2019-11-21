//
//  CCCpy.cpp
//  CCC
//
//  Created by Chris Kiefer on 06/10/2019.
//

/*
https://github.com/TNG/boost-python-examples
*/

#include <stdio.h>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "ETC.hpp"
#include "CCC.hpp"


namespace np = boost::python::numpy;
using namespace boost::python;


char const* greet()
{
   return "Cheese!";
}

double ETC(np::ndarray const & ar) {
  int64_t* input_ptr = reinterpret_cast<int64_t*>(ar.get_data());
  ivec inputVec(input_ptr, ar.shape(0), 1, 1);
  return ETC::calc(inputVec);
}

double dynamicCC(np::ndarray const & ar, size_t dx, size_t xpast, size_t step) {
  int64_t* input_ptr = reinterpret_cast<int64_t*>(ar.get_data());
  ivec inputVec(input_ptr, ar.shape(0), 1, 1);
  return CCC::dynamicCC(inputVec, dx, xpast, step);
}

double dynamicCCJoint(np::ndarray const & ar, np::ndarray const & ar2, size_t dx, size_t past, size_t step) {
  int64_t* input_ptr = reinterpret_cast<int64_t*>(ar.get_data());
  ivec inputVec(input_ptr, ar.shape(0), 1, 1);
  int64_t* input_ptr2 = reinterpret_cast<int64_t*>(ar2.get_data());
  ivec inputVec2(input_ptr2, ar.shape(0), 1, 1);
  return CCC::dynamicCCJoint(inputVec, inputVec2, dx, past, step);
}

double CCCausality(np::ndarray const & ar, np::ndarray const & ar2, size_t dx, size_t past, size_t step) {
  int64_t* input_ptr = reinterpret_cast<int64_t*>(ar.get_data());
  ivec inputVec(input_ptr, ar.shape(0), 1, 1);

  int64_t* input_ptr2 = reinterpret_cast<int64_t*>(ar2.get_data());
  ivec inputVec2(input_ptr2, ar2.shape(0), 1, 1);

  return CCC::CCCausality(inputVec, inputVec2, dx, past, step);
}

BOOST_PYTHON_MODULE(cccrt)
{
    np::initialize();

    def("greet", greet);
    def("ETC", ETC);
    def("dynamicCC", dynamicCC);
    def("dynamicCCJoint", dynamicCCJoint);
    def("CCCausality", CCCausality);
}
