
# Contributions

We encourage contributions to the NetworKit source code. If you encounter any questions regarding NetworKit development please contact the [mailing list][list] or open issues in this GitHub repository.

## Coding Style

- Adapt to the existing naming conventions. NetworKit generally uses lowerCamelCase
for (global and local) variables and  UpperCamelCase for classes.
- We use 4 space indentation.
- Provide tests for new code or modified behavior.

## Usage of Git

- Properly rebase your commits on the development branch before posting or updating a PR.
- Force push to feature branches of your fork to update PRs.
- We prefer a linear history. Do not merge from development branches into feature branches/PRs
  for no good reason.

## Unit tests

NetworKit comes with a large suite of unit tests.
The unit tests can only be run from a clone or copy of the repository and not from a pip installation.
In order to run the unit tests, you need to compile them first.
This is done by setting the [CMake] `NETWORKI_BUILD_TESTS` flag to `ON`:

	cmake -DNETWORKIT_BUILD_TESTS=ON ..

Unit tests are implemented using GTest macros such as `TEST_F(CentralityGTest, testBetweennessCentrality)`.
Single tests can be executed with:

	./networkit_tests --gtest_filter=CentralityGTest.testBetweennessCentrality

Additionally, one can specify the level of the logs outputs by adding `--loglevel=<log_level>`;
supported log levels are: `TRACE`, `DEBUG`, `INFO`, `WARN`, `ERROR`, and `FATAL`.

