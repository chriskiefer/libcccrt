cmake_minimum_required (VERSION 3.15)

project (cccrt)




include_directories(../eigen-3.4.0 ../EigenRand ..)

add_executable(cccrtbench main.cpp ../shannonEntropy.hpp ../ETC.hpp ../CCC.hpp ../RPC.hpp ../fractal.hpp maximilian.cpp)
# target_link_libraries(cccrt gtest_main)



set_property(TARGET cccrtbench PROPERTY CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# add_custom_target(test COMMAND echo "test")

