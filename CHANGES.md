
# API Breaking Changes

API-breaking changes are extracted and filtered by using [ABI-compliance checker](https://lvc.github.io/abi-compliance-checker/) together with [ABI-dumper](https://github.com/lvc/abi-dumper) for C++ core-library and the Cython-interface. Native Python-code functions are covered by diffs.


## Changes in NetworKit 7.0

Due to larger refactoring in the C++ core, there exists quite some changes in source code since the last major release. According to ABI-compliance, there is 96,5% source compatibility between 6.0 to 7.0, so most of the API-calls remain stable.


### Function/Class removed

- Centrality
  - `ApproximatePageRank::push()`


### Function/Class declaration changed

- Distance
  - `SSSP::getPredecessors()` now returns `const std::vector<node>&` instead of `std::vector<node>`


### Changes in Python interface

- `runIteratively()` and `runRecursively()` are remove from `StronglyConnectedComponents`. The default option is now to run in iterative mode. If the (now deprecated) recursive algorithm should be used, the constructor can be called accordingly.
