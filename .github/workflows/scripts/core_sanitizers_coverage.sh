#!/bin/bash

python3 -m venv pyenv && . pyenv/bin/activate
pip3 install --upgrade pip
pip3 install setuptools cpp-coveralls

mkdir build && cd "$_"
cmake -GNinja -DNETWORKIT_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug -DNETWORKIT_WITH_SANITIZERS=leak -DNETWORKIT_COVERAGE=ON ..
ninja

ctest -V