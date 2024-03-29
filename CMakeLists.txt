cmake_minimum_required (VERSION 3.15)

project (cccrt)




# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()


#find_package(Armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})

include_directories(eigen-3.4.0 EigenRand)

# target_include_directories(
#   "${eigen_SOURCE_DIR}"
# )

add_executable(cccrt main.cpp shannonEntropy.hpp ETC.hpp CCC.hpp RPC.hpp fractal.hpp)
#target_link_libraries(cccrt ${ARMADILLO_LIBRARIES})
target_link_libraries(cccrt gtest_main)


#add_test(NAME example_test COMMAND example)

set_property(TARGET cccrt PROPERTY CXX_STANDARD 17)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_custom_target(test COMMAND echo "test")

#project (CCCPython)
#add_library(CCCPython CCCpy.cpp shannonEntropy.hpp ETC.hpp)
#target_link_libraries(CCCPython ${ARMADILLO_LIBRARIES})
