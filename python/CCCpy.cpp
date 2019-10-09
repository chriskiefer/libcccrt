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


BOOST_PYTHON_MODULE(cccrt)
{
    np::initialize();

    def("greet", greet);
    def("ETC", ETC);
    // def("nptest", nptest);
}
