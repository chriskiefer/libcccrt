#!/bin/sh
# Cross-compile the Max externals for Windows from Linux/macOS with MinGW-w64.
#
#   Debian/Ubuntu: sudo apt install g++-mingw-w64-x86-64
#   macOS:         brew install mingw-w64
#
# usage: plugins/max/build-mingw.sh [output-dir]
#   MAX_SDK_BASE_PATH  path to max-sdk-base (default: ../max-sdk-base beside the repo; cloned if missing)
#   output-dir         where to put the .mxe64 files (default: build-mingw/externals under plugins/max)
set -e
here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../.." && pwd)"

sdk="${MAX_SDK_BASE_PATH:-$repo/../max-sdk-base}"
if [ ! -f "$sdk/script/max-pretarget.cmake" ]; then
    echo "cloning max-sdk-base into $sdk"
    git clone --depth 1 https://github.com/Cycling74/max-sdk-base.git "$sdk"
fi
sdk="$(cd "$sdk" && pwd)"

out="${1:-$here/build-mingw/externals}"
mkdir -p "$out"
out="$(cd "$out" && pwd)"

for ext in "cccrpc~" "cccrpc"; do
    build="$here/build-mingw/$ext"
    cmake -S "$here/$ext" -B "$build" \
        -DCMAKE_TOOLCHAIN_FILE="$here/toolchain-mingw-w64.cmake" \
        -DMAX_SDK_BASE_PATH="$sdk" \
        -DC74_LIBRARY_OUTPUT_DIRECTORY="$out" \
        -DCMAKE_BUILD_TYPE=Release -Wno-dev > /dev/null
    cmake --build "$build" 2>&1 | grep -E "error|Error|Built target" || true
done
echo "externals in $out:"
ls -la "$out"
