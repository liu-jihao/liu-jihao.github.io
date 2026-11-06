#!/bin/sh
# Builds viprchk against GMP. Set GMP_PREFIX if GMP is not in a standard location
# (for instance GMP_PREFIX=/opt/homebrew on macOS with Homebrew).
set -e
cd "$(dirname "$0")"
G=${GMP_PREFIX:+-I$GMP_PREFIX/include}
L=${GMP_PREFIX:+-L$GMP_PREFIX/lib}
${CXX:-c++} -O2 -DNDEBUG -std=c++14 -I. $G viprchk.cpp gmp_stream_shim.cpp $L -lgmp -o "${1:-viprchk.out}"
