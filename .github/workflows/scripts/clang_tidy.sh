#!/bin/bash
set -e
set -o pipefail

LLVM_VERSION="${LLVM_VERSION:-22}"

sudo rm -rf /usr/local/clang-*
sudo update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-"${LLVM_VERSION}" 9999
sudo update-alternatives --install /usr/bin/clang clang /usr/bin/clang-"${LLVM_VERSION}" 9999
sudo update-alternatives --install /usr/bin/clang-tidy clang-tidy /usr/bin/clang-tidy-"${LLVM_VERSION}" 9999
$CXX --version
clang-tidy --version
cmake --preset clang-tidy -DCMAKE_CXX_COMPILER="$CXX" -DNETWORKIT_CLANG_TIDY_EXECUTABLE=clang-tidy-"${LLVM_VERSION}"
cmake --build --preset clang-tidy
ctest --preset clang-tidy -V
