#!/bin/bash

python3 -m venv pyenv && . pyenv/bin/activate
pip3 install --upgrade pip
pip3 install setuptools cpp-coveralls

mkdir build && cd "$_"
export CPU_COUNT=$(python3 $GITHUB_WORKSPACE/.github/workflows/scripts/get_core_count.py)
cmake -DNETWORKIT_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug -DNETWORKIT_WITH_SANITIZERS=leak -DNETWORKIT_COVERAGE=ON ..
make -j$CPU_COUNT

ctest -V