<p align="center">
  <img width="60%" src="docs/logo/logo_color.png" alt="NetworKit - Lage-scale Network Analysis"><br>
  <a href="https://github.com/networkit/networkit/actions"><img src="https://github.com/networkit/networkit/workflows/build/badge.svg"></a>
  <a href="https://badge.fury.io/py/networkit"><img src="https://badge.fury.io/py/networkit.svg"></a>
  <a href="https://coveralls.io/github/networkit/networkit?branch=master"><img src="https://coveralls.io/repos/github/networkit/networkit/badge.svg?branch=master"></a>
  <a href="https://mybinder.org/v2/gh/networkit/networkit/master?urlpath=lab/tree/notebooks/User-Guide.ipynb"><img src="https://mybinder.org/badge_logo.svg"></a>
</p>

## 
[NetworKit][networkit] is an open-source tool suite for high-performance
network analysis. Its aim is to provide tools for the analysis of large
networks in the size range from thousands to billions of edges. For this
purpose, it implements efficient graph algorithms, many of them parallel to
utilize multicore architectures. These are meant to compute standard measures
of network analysis. NetworKit is focused on scalability and comprehensiveness.
NetworKit is also a testbed for algorithm engineering and
contains novel algorithms from recently published research (see list of publications below).

NetworKit is a Python module. High-performance algorithms are written in C++ and exposed to Python
via the Cython toolchain. Python in turn gives us the ability to work interactively and a
rich environment of tools for data analysis and scientific computing.
Furthermore, NetworKit's core can be built and used as a native library if needed.

## Requirements

You will need the following software to install NetworKit as a python
package:

- A modern C++ compiler, e.g.: [g++] (&gt;= 6.1), [clang++] (&gt;= 3.9) or MSVC (&gt;= 14.13)
- OpenMP for parallelism (usually ships with the compiler)
- Python3 (3.6 or higher is supported)
  - Development libraries for Python3. The package name depends on your distribution. Examples: 
    - Debian/Ubuntu: `apt-get install python3-dev`
    - RHEL/CentOS: `dnf install python3-devel`
    - Windows: Use the official release installer from [www.python.org](https://www.python.org/downloads/windows/)
- [Pip]
- [CMake] version 3.6 or higher (Advised to use system packages if available. Alternative: `pip3 install cmake`)
- Build system: [Make] or [Ninja]
- Cython version 0.29 or higher (e.g., `pip3 install cython`)

## Install

In order to use NetworKit, you can either install it via package managers
or build the Python module from source.

### Install via package manager

While the most recent version is in general available for all package managers, the number of older downloadable versions differ. 

##### pip

    pip3 install [--user] networkit

##### conda (channel conda-forge)

    conda config --add channels conda-forge
    conda install networkit [-c conda-forge]

##### brew

    brew install networkit

##### spack

    spack install py-networkit

### Building the Python module from source

    git clone https://github.com/networkit/networkit networkit
    cd networkit
    python3 setup.py build_ext [-jX]
    pip3 install -e .

The script will call `cmake` and `ninja` (`make` as fallback) to compile
NetworKit as a library, build the extensions and copy it to the top folder. By
default, NetworKit will be built with the amount of available cores in
optimized mode. It is possible the add the option `-jN` the number of threads
used for compilation.

## Usage example

To get an overview and learn about NetworKit's different functions/classes, have a look at our interactive [notebooks-section][notebooks], especially the [Networkit UserGuide]. Note: To view and edit the computed output from the notebooks, it is recommended to use [Jupyter Notebook][jupyter-notebooks]. This requires the prior installation of NetworKit. You should really check that out before start working on your network analysis. 

We also provide a Binder-instance of our notebooks. To access this service, you can either click on the badge at the top or follow this [link][binder]. Disclaimer: Due to rebuilds of the underlying image, it can takes some time until your Binder instance is ready for usage.

If you only want to see in short how NetworKit is used - the following example provides a climpse at that. Here we generate a random hyperbolic graph with 100k nodes and compute its communities with the PLM method:

    >>> import networkit as nk
    >>> g = nk.generators.HyperbolicGenerator(1e5).generate()
    >>> communities = nk.community.detectCommunities(g, inspect=True)
    PLM(balanced,pc,turbo) detected communities in 0.14577102661132812 [s]
    solution properties:
    -------------------  -----------
    # communities        4536
    min community size      1
    max community size   2790
    avg. community size    22.0459
    modularity              0.987243
    -------------------  -----------

## Install the C++ Core only

In case you only want to work with NetworKit's C++ core, you can either install it via package 
managers or build it from source. 

### Install C++ core via package manager

##### conda (channel conda-forge)

    conda config --add channels conda-forge
    conda install libnetworkit [-c conda-forge]

##### brew

    brew install libnetworkit

##### spack

    spack install libnetworkit

### Building the C++ core from source 

We recommend [CMake] and your preferred build system for building the C++ part of NetworKit.

The following description shows how to use [CMake] in order to build the C++ Core only:

First you have to create and change to a build directory: (in this case named `build`)

    mkdir build
    cd build

Then call [CMake] to generate files for the `make` build system, specifying the directory of the root `CMakeLists.txt` file (e.g., `..`). After this `make` is called to start the build process:

    cmake ..
    make -jX

To speed up the compilation with make a multi-core machine, you can append `-jX` where X denotes the number of threads to compile with.

### Use NetworKit as a library

This paragraph explains how to use the NetworKit core C++ library in case it has been built from source. 
For how to use it when installed via package managers, best refer to the official documentation ([brew](https://brew.sh), [conda](https://docs.conda.io), [spack](https://spack.readthedocs.io/en/latest)). 

In order to use the previous compiled networkit library, you need to have it installed, and link
it while compiling your project. Use these instructions to compile and install NetworKit in `/usr/local`:

    cmake ..
    make -jX install

Once NetworKit has been installed, you can use include directives in your
C++\-application as follows:

    #include <networkit/graph/Graph.hpp>

You can compile your source as follows:

    g++ my_file.cpp -lnetworkit


### Unit tests

Building and running NetworKit unit tests is not mandatory.
However, as a developer you might want to write and run unit tests for your
code, or if you experience any issues with NetworKit, you might want to check
if NetworKit runs properly.
The unit tests can only be run from a clone or copy of the repository and not
from a pip installation. In order to run the unit tests, you need to compile
them first. This is done by setting the [CMake] `NETWORKI_BUILD_TESTS` flag to
`ON`:

    cmake -DNETWORKIT_BUILD_TESTS=ON ..

Unit tests are implemented using GTest macros such as `TEST_F(CentralityGTest, testBetweennessCentrality)`.
Single tests can be executed with:

    ./networkit_tests --gtest_filter=CentralityGTest.testBetweennessCentrality

Additionally, one can specify the level of the logs outputs by adding `--loglevel <log_level>`;
supported log levels are: `TRACE`, `DEBUG`, `INFO`, `WARN`, `ERROR`, and `FATAL`.


### Compiling with address/leak sanitizers

Sanitizers are great tools to debug your code. NetworKit provides additional [Cmake] flags
to enable address, leak, and undefined behavior sanitizers.
To compile your code with sanitizers, set the [CMake]
`NETWORKIT_WITH_SANITIZERS` to either `address` or `leak`:

    cmake -DNETWORKIT_WITH_SANITIZERS=leak ..

By setting this flag to `address`, your code will be compiled with the `address` and the `undefined` sanitizers.
Setting it to `leak` also adds the `leak` sanitizer.


## Documentation

The most recent version of the [documentation can be found online](https://networkit.github.io/dev-docs/index.html).

## Contact

For questions regarding NetworKit, have a look at our [issues-section][issues] and see if there is already an open discussion. If not feel free to open a new issue.
To stay updated about this project, subscribe to our [mailing list][list].

## Contributions

We encourage contributions to the NetworKit source code. See the [development guide][devguide] for instructions. For support please contact the [mailing list][list].

## Credits

List of contributors can be found on the [NetworKit website credits page](https://networkit.github.io/credits.html).

## External Code

The program source includes:
- the *[TLX][tlx]* library
- the *[TTMath]* bignum library

[mitlicense]: http://opensource.org/licenses/MIT
[ttmath]: http://www.ttmath.org/
[tlx]: https://github.com/tlx/tlx/

## License
The source code of this program is released under the [MIT License][mitlicense].  We ask you to cite us if you use this code in your project (c.f. the publications section below and especially the [technical report](https://arxiv.org/abs/1403.3005)). Feedback is also welcome.

## Publications
The [NetworKit publications page][nwkpubs] lists the publications on NetworKit as a toolkit, on algorithms available
in NetworKit, and simply using NetworKit. We ask you to cite the appropriate ones if you found NetworKit useful for your own research.

[nwkpubs]: https://networkit.github.io/publications.html
[list]: https://sympa.cms.hu-berlin.de/sympa/subscribe/networkit
[networkit]: https://networkit.github.io/
[IPython]: https://ipython.readthedocs.io/en/stable/
[NetworKit UserGuide]: https://github.com/networkit/networkit/blob/master/notebooks/User-Guide.ipynb
[notebooks]: https://github.com/networkit/networkit/blob/master/notebooks/
[g++]: https://gcc.gnu.org
[clang++]: https://clang.llvm.org/
[Pip]: https://pypi.python.org/pypi/pip
[CMake]: https://cmake.org/
[Make]: https://www.gnu.org/software/make/
[Ninja]: https://ninja-build.org/
[devguide]: https://networkit.github.io/dev-docs/DevGuide.html#devGuide
[issues]: https://github.com/networkit/networkit/issues
[jupyter-notebooks]: https://jupyter.org/install.html
[binder]: https://mybinder.org/v2/gh/networkit/networkit/master?urlpath=lab/tree/notebooks
