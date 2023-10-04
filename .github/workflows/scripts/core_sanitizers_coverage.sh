#!/bin/bash

python3 -m venv pyenv && . pyenv/bin/activate
pip3 install --upgrade pip
pip3 install coveralls cython gcovr matplotlib requests setuptools tabulate

mkdir core_build && cd "$_"
export CPU_COUNT=$(python3 $GITHUB_WORKSPACE/.github/workflows/scripts/get_core_count.py)
cmake -DNETWORKIT_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug -DNETWORKIT_WITH_SANITIZERS=leak -DNETWORKIT_COVERAGE=ON ..
make -j$CPU_COUNT

ctest -V

cd ..
export CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH:+$CMAKE_LIBRARY_PATH:}$(pwd)/core_build
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}$(pwd)/core_build

NETWORKIT_PARALLEL_JOBS=$CPU_COUNT python3 setup.py build_ext --inplace
NETWORKIT_PARALLEL_JOBS=$CPU_COUNT pip3 install -e .