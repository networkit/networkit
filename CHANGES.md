
# API Breaking Changes

API-breaking changes are extracted and filtered by using [ABI-compliance checker](https://lvc.github.io/abi-compliance-checker/) together with [ABI-dumper](https://github.com/lvc/abi-dumper) for C++ core-library and the Cython-interface. Native Python-code functions are covered by diffs.

## Changes in NetworKit 8.0

### Function/Class removed

- Aux
  - `Parrallelism::enableNestedParallelism()` 

- Centrality
  - Remark: The implementation of `TopHarmonicCloseness` was reworked entirely. In the new version all formerly `protected` functions are removed. The `public` API has not changed.
  - `TopHarmonicCloseness::BFSbound()`
  - `TopHarmonicCloseness::BFScut()`
  - `TopHarmonicCloseness::computeReachableNodesDirected()`
  - `TopHarmonicCloseness::computeReachableNodesUndirected()`
  - `TopHarmonicCloseness::init()`
  - `Centrality::scores( bool moveOut )`

- Components
  - Remark: The recursive version of `StronglyConnectedComponents` is removed. To call the algorithm use `StronglyConnectedComponents::run()`.
  - `StronglyConnectedComponents::runIteratively()`
  - `StronglyConnectedComponents::runRecursively()`
  - `StronglyConnectedComponents::StronglyConnectedComponents( Graph const& G, bool iterativeAlgo )`

- In class `Graph` quite some functions are removed with release of 8.0. They are now mostly available via `GraphTools`:
  - `Graph::append()`
  - `Graph::copyNodes()`
  - `Graph::density()`
  - `Graph::edges()`
  - `Graph::maxDegree()`
  - `Graph::maxDegreeIn()`
  - `Graph::maxWeightedDegree()`
  - `Graph::maxWeightedDegreeIn()`
  - `Graph::merge()`
  - `Graph::neighbors()`
  - `Graph::nodes()`
  - `Graph::randomEdge()`
  - `Graph::randomEdges()`
  - `Graph::randomNeighbor()`
  - `Graph::randomNode()`
  - `Graph::removeEdgesFromIsolatedSet()`
  - `Graph::size()`
  - `Graph::subgraphFromNodes()`
  - `Graph::toString()`
  - `Graph::toUndirected()`
  - `Graph::toUnweighted()`
  - `Graph::transpose()`
  - `Graph::typ()`
  - `Graph::volume()`


### Function/Class declaration changed

- Distances
  - `SSSP::getDistances( bool moveOut )` to `SSSP::getDistances()`
  - `SSSP::getNodesSortedByDistance( bool moveOut )` to `SSSP::getNodesSortedByDistance()`


### Function/Class now deprecated

- Randomization
  - `GlobalTradeSequence.hpp` is now deprecated

### Changes in Python interface

- components
  - `StronglyConnectedComponents` no longer supports recursive solving strategy.
- engineering
  - `enableNestedParallelism()` is no longer available
- graph
  - The following functions are removed from the `graph`-module. Most of them are available now from the `graphtools`-module:
  `append()`, `copyNodes()`, `density()`, `edges()`, `inNeighbors`, `maxDegree()`, `maxDegreeIn()`, `maxWeightedDegree()`, `maxWeightedDegreeIn()`, `merge()`, `neighbors()`, `randomEdge()`, `randomEdges()`, `randomNeighbor()`, `randomNode()`, `removeEdgesFromIsolatedSet()`, `size()`, `subgraphFromNodes()`, `toString()`, `toUndirected()`, `toUnweighted()`, `transpose()`, `typ()`, `volume()`

## Changes in NetworKit 7.0

### Function/Class removed

- Centrality
  - `ApproximatePageRank::push()`


### Function/Class declaration changed

- Distance
  - `SSSP::getPredecessors()` now returns `const std::vector<node>&` instead of `std::vector<node>`


### Changes in Python interface

- `runIteratively()` and `runRecursively()` are removed from `StronglyConnectedComponents`. The default option is now to run in iterative mode. If the (now deprecated) recursive algorithm should be used, the constructor can be called accordingly.
