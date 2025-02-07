#!/bin/bash
set -e
set -o pipefail

$CXX --version
python3 --version
cmake --version

# Note: setuptools<69.0.0 is needed for editable installs. Normal installs work with more recent releases of setuptools
python3 -m venv pyenv && . pyenv/bin/activate
pip3 install --upgrade pip
pip3 install cython "numpy" ipython jupyter setuptools

# Build tlx
cd tlx
mkdir build && cd "$_"
cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$TLX_PATH ..
ninja
ninja install
cd ../..


# Build libnetworkit
mkdir core_build && cd "$_"
cmake -GNinja -DNETWORKIT_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release -DNETWORKIT_CXX_STANDARD=$CXX_STANDARD -DNETWORKIT_WARNINGS=ON -DNETWORKIT_WARNINGS_AS_ERRORS=ON -DNETWORKIT_EXT_TLX=$TLX_PATH -DNETWORKIT_NATIVE=$NATIVE ..
ninja
UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ASAN_OPTIONS=abort_on_error=1 ctest -V
cd ..
export CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH:+$CMAKE_LIBRARY_PATH:}$(pwd)/core_build
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}$(pwd)/core_build
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH:+$DYLD_LIBRARY_PATH:}$(pwd)/core_build
export CPU_COUNT=$(python3 $GITHUB_WORKSPACE/.github/workflows/scripts/get_core_count.py)

# Build networkit
NETWORKIT_PARALLEL_JOBS=$CPU_COUNT python3 ./setup.py build_ext --inplace --networkit-external-core --external-tlx=$TLX_PATH
NETWORKIT_PARALLEL_JOBS=$CPU_COUNT pip3 install -e .

# Test Python module
python3 -c 'import networkit'

pip3 install -r requirements.txt

python3 -m unittest discover -v networkit/test/

# Installing ipywidgets and jupyterlab-widgets (+ pinned versions) is temp. necessary due to incompatibility with Jupyterlab/plotly and ipywidgets>8.X
# See: https://github.com/plotly/plotly.py/issues/3686
pip3 install ipycytoscape plotly seaborn ipywidgets==7.7.1 jupyterlab-widgets==1.1.1 anywidget

python3 notebooks/test_notebooks.py 'notebooks/'
