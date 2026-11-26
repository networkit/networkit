#!/bin/bash
set -e
set -o pipefail

$CXX --version

cmake --preset debug -DNETWORKIT_MONOLITH=$MONOLITH -DNETWORKIT_CXX_STANDARD=$CXX_STANDARD
cmake --build --preset debug

ctest --preset debug -V
