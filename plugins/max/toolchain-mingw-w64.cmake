# CMake toolchain for cross-compiling Max externals from Linux/macOS with MinGW-w64.
#   Debian/Ubuntu: apt install g++-mingw-w64-x86-64
#   macOS:         brew install mingw-w64
# usage: cmake -S cccrpc~ -B build-win -DCMAKE_TOOLCHAIN_FILE=../toolchain-mingw-w64.cmake -DMAX_SDK_BASE_PATH=...
set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_SYSTEM_PROCESSOR AMD64)

set(TOOLCHAIN_PREFIX x86_64-w64-mingw32)
set(CMAKE_C_COMPILER   ${TOOLCHAIN_PREFIX}-gcc)
set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)
set(CMAKE_RC_COMPILER  ${TOOLCHAIN_PREFIX}-windres)

set(CMAKE_FIND_ROOT_PATH /usr/${TOOLCHAIN_PREFIX})
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

# the Max SDK headers are written for MSVC; GCC just needs to ignore its #pragma warning lines
set(CMAKE_CXX_FLAGS_INIT "-Wno-unknown-pragmas")
set(CMAKE_C_FLAGS_INIT "-Wno-unknown-pragmas")
