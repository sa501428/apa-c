#!/bin/bash

set -euo pipefail

ZSTD_CFLAGS="$(pkg-config --cflags libzstd 2>/dev/null || true)"
ZSTD_LIBS="$(pkg-config --libs libzstd 2>/dev/null || true)"

g++ -std=c++17 -O3 -Wall -Wextra \
  ${ZSTD_CFLAGS} \
  -I../straw-master/C++ \
  -o apa_hic \
  main.cpp apa_hic.cpp ../straw-master/C++/straw.cpp \
  -lz -pthread \
  ${ZSTD_LIBS}
