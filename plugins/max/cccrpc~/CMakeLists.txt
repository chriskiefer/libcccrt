cmake_minimum_required(VERSION 3.19)

# Path to a checkout of https://github.com/Cycling74/max-sdk-base
# (or the max-sdk-base folder inside a max-sdk checkout: max-sdk/source/max-sdk-base)
if (NOT DEFINED MAX_SDK_BASE_PATH)
	set(MAX_SDK_BASE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../../max-sdk-base")
endif ()
if (NOT EXISTS "${MAX_SDK_BASE_PATH}/script/max-pretarget.cmake")
	message(FATAL_ERROR "max-sdk-base not found at ${MAX_SDK_BASE_PATH}; set -DMAX_SDK_BASE_PATH=/path/to/max-sdk-base")
endif ()

# where the built external goes (must be set before max-pretarget)
set(C74_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../externals")

include(${MAX_SDK_BASE_PATH}/script/max-pretarget.cmake)

include_directories(
	"${MAX_SDK_INCLUDES}"
	"${MAX_SDK_MSP_INCLUDES}"
	"${MAX_SDK_JIT_INCLUDES}"
	"${CMAKE_CURRENT_SOURCE_DIR}/../../../core"
)

add_library(${PROJECT_NAME} MODULE "${CMAKE_CURRENT_SOURCE_DIR}/cccrpc~.cpp")

include(${MAX_SDK_BASE_PATH}/script/max-posttarget.cmake)
