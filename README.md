<div id="top"></div>

<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://networkit.github.io/">
    <img width="60%" src="docs/logo/logo_color.png" alt="NetworKit - Large-scale Network Analysis">
  </a>

  <h3 align="center">NetworKit</h3>

  <p align="center">
    High-performance tools for large-scale network analysis
    <br />
    <a href="https://networkit.github.io/dev-docs/index.html"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://mybinder.org/v2/gh/networkit/networkit/master?urlpath=lab/tree/notebooks/User-Guide.ipynb">Try Demo</a>
    ·
    <a href="https://github.com/networkit/networkit/issues">Report Bug</a>
    ·
    <a href="https://github.com/networkit/networkit/issues">Request Feature</a>
  </p>

  <p align="center">
    <a href="https://github.com/networkit/networkit/actions"><img src="https://github.com/networkit/networkit/workflows/build/badge.svg"></a>
    <a href="https://badge.fury.io/py/networkit"><img src="https://badge.fury.io/py/networkit.svg"></a>
    <a href="https://coveralls.io/github/networkit/networkit?branch=master"><img src="https://coveralls.io/repos/github/networkit/networkit/badge.svg?branch=master"></a>
    <a href="https://mybinder.org/v2/gh/networkit/networkit/master?urlpath=lab/tree/notebooks/User-Guide.ipynb"><img src="https://mybinder.org/badge_logo.svg"></a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#key-features">Key Features</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#building-from-source">Building from Source</a></li>
    <li><a href="#documentation">Documentation</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#publications">Publications</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

[NetworKit][networkit] is an open-source toolkit for high-performance network analysis, designed to handle large networks ranging from thousands to billions of edges. Built with efficiency and scalability at its core, NetworKit implements parallel graph algorithms that leverage multicore architectures to compute standard measures of network analysis.

As both a production tool and a research testbed for algorithm engineering, NetworKit includes novel algorithms from recent publications alongside battle-tested implementations. The toolkit is available as a Python module with high-performance C++ algorithms exposed through Cython, combining Python's interactivity and rich ecosystem with C++'s computational efficiency.

<p align="right">(<a href="#top">back to top</a>)</p>

### Key Features

- **Scalable**: Analyze networks with billions of edges
- **Fast**: Parallel algorithms utilizing multicore architectures
- **Comprehensive**: Wide range of network analysis algorithms
- **Interactive**: Python interface with Jupyter notebook support
- **Flexible**: Available as Python module or standalone C++ library
- **Research-ready**: Includes state-of-the-art algorithms from recent publications

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

### Installation

#### Python Module

**For most users**, NetworKit can be installed directly via package managers with no additional requirements other than Python 3.9+.

| Package Manager | Command |
|-----------------|---------|
| **pip** | `pip install networkit` |
| **conda** | `conda install -c conda-forge networkit` |
| **brew** | `brew install networkit` |
| **spack** | `spack install py-networkit` |

#### C++ Core Library Only

If you only need the C++ core without Python bindings:

| Package Manager | Command |
|-----------------|---------|
| **conda** | `conda install -c conda-forge libnetworkit` |
| **brew** | `brew install libnetworkit` |
| **spack** | `spack install libnetworkit` |

More platform-specific installation instructions can be found in our [getting started guide](https://networkit.github.io/get_started.html).

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

Here's a quick example showing how to generate a random hyperbolic graph with 100k nodes and detect communities:

```python
from networkit.generators import HyperbolicGenerator
from networkit.community import detectCommunities

# Generate a random hyperbolic graph
g = (
    HyperbolicGenerator(1e5)
    .generate()
)

# Detect communities
detectCommunities(g, inspect=True)
```

Output:
```
PLM(balanced,pc,turbo) detected communities in 0.14577102661132812 [s]
solution properties:
-------------------  -----------
# communities        4536
min community size      1
max community size   2790
avg. community size    22.0459
modularity              0.987243
-------------------  -----------
```

### More Examples

Compute PageRank to rank nodes by importance:

```python
from networkit.centrality import PageRank

pr = (
    PageRank(g)
    .run()
)
top_nodes = pr.ranking()[:10]
```

Analyze graph structure with connected components:

```python
from networkit.components import ConnectedComponents

cc = (
    ConnectedComponents(g)
    .run()
)
print(f"Components: {cc.numberOfComponents()}")
print(f"Largest: {max(cc.getComponentSizes().values())}")
```


For comprehensive examples and tutorials, explore our [interactive notebooks][notebooks], especially the [NetworKit User Guide][NetworKit UserGuide]. You can try NetworKit directly in your browser using our [Binder instance][binder].

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- BUILDING FROM SOURCE -->
## Building from Source

### Prerequisites

Building from source requires:

- **C++ Compiler**: [g++] (&gt;= 10.0), [clang++] (&gt;= 11.0), or MSVC (&gt;= 14.30)
- **OpenMP**: For parallelism (usually included with compiler)
- **Python**: 3.9 or higher with development libraries
  - Debian/Ubuntu: `apt-get install python3-dev`
  - RHEL/CentOS: `dnf install python3-devel`
  - Windows: [Official installer](https://www.python.org/downloads/windows/)
- **[CMake]**: Version 3.6 or higher
- **Build System**: [Make] or [Ninja]

### Python Module

```bash
git clone https://github.com/networkit/networkit networkit
cd networkit
pip install cython numpy setuptools wheel
python setup.py build_ext [-jX]
pip install -e .
```

The `-jX` option specifies the number of threads for compilation (e.g., `-j4` for 4 threads). If omitted, it uses all available CPU cores.

### C++ Core Library

```bash
mkdir build && cd build
cmake ..
make -jX
sudo make install
```

#### Using NetworKit in Your C++ Project

After installation, include NetworKit headers:

```cpp
#include <networkit/graph/Graph.hpp>
```

Compile your project:

```bash
g++ my_file.cpp -lnetworkit
```

#### Running Unit Tests

To build and run tests:

```bash
cmake -DNETWORKIT_BUILD_TESTS=ON ..
make
./networkit_tests --gtest_filter=CentralityGTest.testBetweennessCentrality
```

#### Building with Sanitizers

For debugging with address/leak sanitizers:

```bash
cmake -DNETWORKIT_WITH_SANITIZERS=leak ..
```

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- DOCUMENTATION -->
## Documentation

The complete documentation is available online at [networkit.github.io](https://networkit.github.io/dev-docs/index.html).

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

We welcome contributions to NetworKit! Whether you're fixing bugs, adding features, or improving documentation, your help makes NetworKit better for everyone.

1. Check our [development guide][devguide] for instructions
2. Browse [open issues][issues] or open a new one
3. Fork the repository and create your feature branch
4. Submit a pull request

For support, join our [mailing list][list].

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- LICENSE -->
## License

Distributed under the [MIT License][mitlicense]. We ask that you cite us if you use NetworKit in your research (see our [technical report](https://arxiv.org/abs/1403.3005) and [publications page][nwkpubs]).

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

- **Issues**: Check our [issues section][issues] for existing discussions or open a new issue
- **Mailing List**: [Subscribe here][list] to stay updated

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- PUBLICATIONS -->
## Publications

NetworKit has been used in numerous research projects. Visit our [publications page][nwkpubs] for a complete list of papers about NetworKit, algorithms implemented in NetworKit, and research using NetworKit.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CREDITS -->
## Credits

NetworKit is developed by a dedicated team of researchers and contributors. View the full list of contributors on our [credits page](https://networkit.github.io/credits.html).

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
[mitlicense]: http://opensource.org/licenses/MIT
[ttmath]: http://www.ttmath.org/
[tlx]: https://github.com/tlx/tlx/
[nwkpubs]: https://networkit.github.io/publications.html
[list]: https://sympa.cms.hu-berlin.de/sympa/subscribe/networkit
[networkit]: https://networkit.github.io/
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
[binder]: https://mybinder.org/v2/gh/networkit/networkit/master?urlpath=lab/tree/notebooks
