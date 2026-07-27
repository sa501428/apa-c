#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
STRAW_CPP_DIR="${STRAW_CPP_DIR:-"${SCRIPT_DIR}/../../straw/C++"}"

if [[ ! -f "${STRAW_CPP_DIR}/straw.cpp" || ! -f "${STRAW_CPP_DIR}/straw.h" ]]; then
  echo "C++ straw was not found at ${STRAW_CPP_DIR}" >&2
  echo "Set STRAW_CPP_DIR to the directory containing straw.cpp and straw.h." >&2
  exit 1
fi

ZSTD_CFLAGS="$(pkg-config --cflags libzstd 2>/dev/null || true)"
ZSTD_LIBS="$(pkg-config --libs libzstd 2>/dev/null || true)"
CURL_CFLAGS="$(pkg-config --cflags libcurl 2>/dev/null || true)"
CURL_LIBS="$(pkg-config --libs libcurl 2>/dev/null || true)"

[[ -n "${ZSTD_LIBS}" ]] || ZSTD_LIBS="-lzstd"
[[ -n "${CURL_LIBS}" ]] || CURL_LIBS="-lcurl"

g++ -std=c++17 -O3 -DNDEBUG -march=native -mtune=native -flto -Wall -Wextra \
  ${ZSTD_CFLAGS} ${CURL_CFLAGS} \
  -I"${STRAW_CPP_DIR}" \
  -o "${SCRIPT_DIR}/apa_hic" \
  "${SCRIPT_DIR}/main.cpp" "${SCRIPT_DIR}/apa_hic.cpp" "${STRAW_CPP_DIR}/straw.cpp" \
  -lz -pthread \
  ${ZSTD_LIBS} ${CURL_LIBS}
