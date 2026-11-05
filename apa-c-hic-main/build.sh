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
HTS_CFLAGS="$(pkg-config --cflags htslib 2>/dev/null || true)"
HTS_LIBS="$(pkg-config --libs htslib 2>/dev/null || true)"
HTS_DEFINE=""
STRAW_EXTRA_SOURCES=()

[[ -n "${ZSTD_LIBS}" ]] || ZSTD_LIBS="-lzstd"
[[ -n "${CURL_LIBS}" ]] || CURL_LIBS="-lcurl"

if [[ -n "${HTS_LIBS}" ]]; then
  HTS_DEFINE="-DHAVE_HTSLIB"
else
  echo "htslib was not found; building without BAM input support" >&2
fi

# Newer straw checkouts split v10/HBS support into separate translation units.
[[ ! -f "${STRAW_CPP_DIR}/straw_v10.cpp" ]] || STRAW_EXTRA_SOURCES+=("${STRAW_CPP_DIR}/straw_v10.cpp")
[[ ! -f "${STRAW_CPP_DIR}/hbs.cpp" ]] || STRAW_EXTRA_SOURCES+=("${STRAW_CPP_DIR}/hbs.cpp")

g++ -std=c++17 -O3 -DNDEBUG -march=native -mtune=native -flto -Wall -Wextra \
  ${ZSTD_CFLAGS} ${CURL_CFLAGS} ${HTS_CFLAGS} ${HTS_DEFINE} \
  -I"${STRAW_CPP_DIR}" \
  -o "${SCRIPT_DIR}/apa_hic" \
  "${SCRIPT_DIR}/main.cpp" "${SCRIPT_DIR}/apa_hic.cpp" "${STRAW_CPP_DIR}/straw.cpp" \
  "${STRAW_EXTRA_SOURCES[@]}" \
  -lz -pthread \
  ${ZSTD_LIBS} ${CURL_LIBS} ${HTS_LIBS}
