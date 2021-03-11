#!/bin/bash
set -e
set -o pipefail

sudo rm -rf /usr/local/clang-*
sudo update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-11 9999
sudo update-alternatives --install /usr/bin/clang clang /usr/bin/clang-11 9999
sudo update-alternatives --install /usr/bin/clang-tidy clang-tidy /usr/bin/clang-tidy-11 9999
$CXX --version
clang-tidy --version
mkdir debug_test && cd "$_"
cmake -GNinja -DNETWORKIT_BUILD_TESTS=ON -DNETWORKIT_CLANG_TIDY=ON -DCMAKE_BUILD_TYPE=Debug ..
ninja
ctest -V
