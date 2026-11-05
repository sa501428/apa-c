#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CXX="${CXX:-c++}"

# Use a source checkout with a built static library, or an installed libBigWig.
BIGWIG_CFLAGS=()
BIGWIG_LIBS=(-lBigWig)
if [[ -n "${LIBBIGWIG_DIR:-}" ]]; then
  if [[ ! -f "${LIBBIGWIG_DIR}/bigWig.h" || ! -f "${LIBBIGWIG_DIR}/libBigWig.a" ]]; then
    echo "LIBBIGWIG_DIR must contain bigWig.h and libBigWig.a (run make lib-static there)." >&2
    exit 1
  fi
  BIGWIG_CFLAGS=(-I"${LIBBIGWIG_DIR}")
  BIGWIG_LIBS=("${LIBBIGWIG_DIR}/libBigWig.a")
elif command -v pkg-config >/dev/null && pkg-config --exists libBigWig; then
  read -r -a BIGWIG_CFLAGS <<< "$(pkg-config --cflags libBigWig)"
  read -r -a BIGWIG_LIBS <<< "$(pkg-config --libs libBigWig)"
fi

CURL_CFLAGS=()
CURL_LIBS=(-lcurl)
if command -v pkg-config >/dev/null && pkg-config --exists libcurl; then
  read -r -a CURL_CFLAGS <<< "$(pkg-config --cflags libcurl)"
  read -r -a CURL_LIBS <<< "$(pkg-config --libs libcurl)"
fi

"${CXX}" -O3 -std=c++17 -Wall -Wextra \
  ${BIGWIG_CFLAGS[@]+"${BIGWIG_CFLAGS[@]}"} ${CURL_CFLAGS[@]+"${CURL_CFLAGS[@]}"} \
  "${SCRIPT_DIR}/ata_bigwig.cpp" -o "${SCRIPT_DIR}/ata_bigwig" \
  "${BIGWIG_LIBS[@]}" "${CURL_LIBS[@]}" -lz -lm
