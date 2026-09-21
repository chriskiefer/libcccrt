#!/bin/sh
# Cross-compile the core tests for Cortex-M33 (RP2350) to confirm the core
# headers build with an embedded toolchain. Needs arm-none-eabi-g++ on PATH.
set -e
cd "$(dirname "$0")"
OUT="${1:-core_test_m33.elf}"
arm-none-eabi-g++ -std=c++17 -O2 -Wall -Wextra \
    -mcpu=cortex-m33 -mthumb -mfloat-abi=hard -mfpu=fpv5-sp-d16 \
    -fno-exceptions -fno-rtti --specs=nosys.specs --specs=nano.specs \
    core_test.cpp -o "$OUT" 2>&1 | grep -v "is not implemented and will always fail" | grep -v "in function \`_" || true
arm-none-eabi-size "$OUT"
