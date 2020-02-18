#  Contribution Guidelines and Policies

## Contribution Guidelines

We encourage contributions to the NetworKit source code.
If you encounter any questions regarding NetworKit development open an issue on GitHub or contact the [mailing list][list].

### General Guidelines

- NetworKit comes with a large collection of unit tests.
  Please provide tests for new code and for behavior that is modified.
  See below for instructions on building and running the tests.
- We generally welcome Jupyter notebooks that demonstrate how new algorithms are used. If you add a new feature,
  consider adding a notebook to the `notebooks/` directory (or add your feature to an existing notebook).

### Coding Style

- We use 4 space indentation.
- Adapt to the existing naming conventions. NetworKit generally uses lowerCamelCase
  for (global and local) variables and UpperCamelCase for classes.
- We have a `.clang-format` file in the root directory for C++ code styling. To use it, call `clang-format -style=file <code-file>`. Note: We do not require all code to follow the coding style specified. However, it should be taken as a strong suggestion. If you use clang-format, add an extra line with `// networkit-format` before the includes.     

###  Usage of Git

We prefer a linear Git history. In particular:
- Properly rebase your commits on the development branch before posting or updating a PR.
- Do not merge from development branches into feature branches/PRs for no good reason.
- Force push to the feature branches of your fork to update PRs.

###  Building and Running Unit tests

NetworKit comes with a large suite of unit tests. We use the Google Test framework to define and run our tests.
The unit tests can only be run from a clone or copy of the repository and not from a pip installation. In order to run the unit tests, you need to compile them first.

This is done by setting the [CMake] `NETWORKI_BUILD_TESTS` flag to `ON`:

    cmake -DNETWORKIT_BUILD_TESTS=ON ..

Single tests can be executed with:

    ./networkit_tests --gtest_filter=CentralityGTest.testBetweennessCentrality

Additionally, one can specify the level of the logs outputs by adding `--loglevel=<log_level>`;
supported log levels are: `TRACE`, `DEBUG`, `INFO`, `WARN`, `ERROR`, and `FATAL`. For more details, see the
output of `./networkit_tests -h`.

## Policies and Maintainance

### Supported Compilers

We aim to support all releases of GCC and Clang for 5 years. In effect, this means that we can only use C++ features
that have been available in a released GCC or Clang compiler for more than 5 years.

### Deprecation and Removal of APIs

Features need to be deprecated for two releases before they are removed.

*Breaking API changes* should be kept to a minimum. If breaking API changes are done, they should be
clearly documented in CHANGES.md (including guidelines on migrating existing code).
Breaking ABI changes are allowed on every major release. NetworKit does not guarantee a stable ABI
across major releases
(i.e., after upgrading the NetworKit library, you need to recompile programs linking against the library).
Minor and patch releases guarantee ABI stability.

### PR Reviews and Merging

- PRs require at least one review (by a maintainer who is not also the author) before they can be merged.
- PRs should be kept open for at least 7 days *after the last substantial change*.
  This allows others to take a look at the latest version of the PR before it is merged.
  Critical bugfixes are exempt from this policy.

[list]: https://sympa.cms.hu-berlin.de/sympa/subscribe/networkit
