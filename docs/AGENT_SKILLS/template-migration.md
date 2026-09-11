---

name: networkit-template-migration
description: Use when migrating a NetworKit C++ algorithm, data structure, helper class, reader, or writer from fixed Graph/node/edgeweight types to the templated graph architecture. Preserve the existing public API, use GraphT/NodeT/EdgeWeightT consistently, generalize tests with typed tests, move template implementations out of source .cpp files, delete obsolete .cpp files and their CMake entries, and avoid performance regressions.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# NetworKit Template Migration

You are an expert C++ contributor working on NetworKit's incremental migration from fixed graph types to a templated graph library.

The goal is not merely to make existing code compile with `AdjListGraph<NodeT, EdgeWeightT>`. The migration must preserve existing behavior and source compatibility while making graph-dependent types genuinely generic.

## Core principles

1. **Preserve the existing API.** Existing user-facing class names should continue to work with the traditional `Graph` configuration.

2. **Introduce a genuinely generic implementation.** Replace hard-coded graph-dependent types with types derived from the template parameters.

3. **Do not change semantics accidentally.** Templating is primarily a type-generalization change. Avoid unrelated behavioral changes, numerical changes, cleanup, or API redesign unless required for correctness.

4. **Use semantic types rather than mechanical replacement.** Do not blindly replace every `node`, `edgeweight`, `index`, or `count`. Determine what each value represents.

5. **Generalize the tests together with the implementation.** A templated implementation without tests covering multiple type combinations is incomplete.

6. **Remove obsolete source translation units.** When an implementation is moved from a `.cpp` file into a template header or `Impl.hpp`, delete the original `.cpp` file and remove it from the corresponding `CMakeLists.txt`. Do not leave an empty, forwarding, or compatibility `.cpp` behind.

7. **Keep the default specialization efficient.** Avoid additional indirection, unnecessary conversions, allocations, or abstractions in hot code. For performance-sensitive data structures or algorithms, compare the default specialization against the previous implementation.

8. **Keep the PR focused.** Rebase onto required prerequisite changes before implementing or reviewing the migration. Do not reintroduce outdated changes from dependent branches.

## Public API and naming

When an existing class must become templated, introduce a new templated implementation name and preserve the old name as an alias.

For algorithms, a typical pattern is:

```cpp
template <class GraphT>
class GenericFoo final : public Algorithm {
public:
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    explicit GenericFoo(const GraphT &graph);

    void run() override;

private:
    const GraphT *graph;
};

using Foo = GenericFoo<Graph>;
```

If the old class name is an abbreviation and there is a natural descriptive name, a descriptive templated name may be preferable:

```cpp
template <class GraphT>
class BreadthFirstSearch : public SingleSourceShortestPaths<GraphT> {
    // ...
};

using BFS = BreadthFirstSearch<Graph>;
```

Follow the naming convention already used by nearby migrated classes. Do not rename the legacy API exposed to existing users.

For generic data structures, template the relevant value types directly:

```cpp
template <SomeKeyConcept KeyType, SomeValueConcept ValueType>
class DataStructure {
    // ...
};

using OldDataStructure = DataStructure<OldKeyType, OldValueType>;
```

When template parameters have restrictions, express them using existing NetworKit concepts or appropriate C++20 concepts instead of relying on obscure template errors.

## Type migration

Start by auditing all occurrences of:

```text
Graph
node
edgeweight
none
index
count
edgeid
numeric literals used as weights/distances
containers storing graph-dependent values
lambda parameters passed to graph iteration functions
```

Classify each occurrence by meaning.

Use:

```cpp
using NodeT = typename GraphT::NodeT;
using EdgeWeightT = typename GraphT::EdgeWeightT;
```

for graph-dependent node identifiers and edge weights.

Examples:

```cpp
const GraphT *graph;

std::vector<NodeT> roots;

graph->forNodes([&](NodeT u) {
    graph->forNeighborsOf(u, [&](NodeT v) {
        // ...
    });
});
```

Use the templated node sentinel:

```cpp
NullNodeId<NodeT>
```

instead of a fixed node sentinel such as `none` when the value represents a node ID.

Use the canonical edge-ID sentinel, such as:

```cpp
nullEdgeId
```

when the value is an edge ID and edge IDs themselves are not part of the graph template parameters.

Do **not** convert every `count` or `index` to `NodeT`.

Keep `count` when the value semantically represents a count and `index` when it represents a generic array/indexing quantity independent of the graph's node-ID representation.

Conversely, variables that actually represent node identifiers, node limits passed through a templated graph API, or node values used to construct test graphs should normally use `NodeT`.

### Distinguish graph types from algorithm result types

Do not assume that every numerical result should become `EdgeWeightT`.

An algorithm may intentionally distinguish:

```cpp
using EdgeWeightT = typename GraphT::EdgeWeightT;
using DistanceT = edgeweight;
```

for example when the historic API promises `edgeweight`, when special floating-point values are required, or when the result type has semantics different from the graph's physical edge-storage type.

Preserve the existing API's numerical semantics unless changing them is an explicit goal.

Make such distinctions explicit through aliases rather than scattered casts.

### Numeric constants

Avoid unnecessarily forcing calculations through `double`.

Prefer type-aware constants:

```cpp
constexpr EdgeWeightT one{1};
distances[source] = EdgeWeightT{0};
```

when the operation really belongs to `EdgeWeightT`.

### Signed and unsigned conversions

Generic types expose sign-comparison problems that fixed types may have hidden.

When a generic integral value must index a vector:

```cpp
if constexpr (std::is_signed_v<ValueType>) {
    assert(value >= 0);
}

const auto valueIndex = static_cast<index>(value);
assert(valueIndex < values.size());
```

Perform conversions at a small number of explicit boundaries. Do not scatter casts throughout the algorithm merely to silence warnings.

If a generic data structure has logically different types — for example key, stored value, and backing-array index — keep those types distinct.

## Template implementation layout

Template definitions must be visible at the point of instantiation. Therefore, implementations that previously lived in `.cpp` files generally have to move into headers.

Use one of two layouts.

For a compact algorithm whose declaration and implementation remain easy to understand, keep everything in the main header:

```text
Foo.hpp
```

For a larger implementation, use:

```text
Foo.hpp
FooImpl.hpp
```

with declarations and public documentation in `Foo.hpp`, implementations in `FooImpl.hpp`, and:

```cpp
#include <networkit/.../FooImpl.hpp>
```

at the bottom of `Foo.hpp`.

Do not create an `Impl.hpp` merely because a class became templated. Choose the split based on readability and the conventions of nearby migrated components.

### Mandatory removal of the original `.cpp`

If the pre-migration implementation lived in:

```text
networkit/cpp/<module>/Foo.cpp
```

and the implementation has been moved into `Foo.hpp` or `FooImpl.hpp`, then:

```text
networkit/cpp/<module>/Foo.cpp
```

**must be deleted.**

Do not:

* leave an empty `.cpp`;
* leave a `.cpp` that only includes `Foo.hpp`;
* keep explicit instantiations solely to preserve the old file;
* keep a forwarding `.cpp` for compatibility;
* leave the file in place because the build still succeeds.

The compatibility mechanism is the legacy type alias in the public header, not a legacy `.cpp` file.

### Mandatory CMake cleanup

Whenever such a `.cpp` file is removed, find the corresponding module `CMakeLists.txt` and remove the source entry in the same change.

For example, convert:

```cmake
networkit_add_module(distance
    Dijkstra.cpp
    FloydWarshall.cpp
    GraphDistance.cpp
)
```

to:

```cmake
networkit_add_module(distance
    Dijkstra.cpp
    GraphDistance.cpp
)
```

The migration is incomplete if:

* the original implementation `.cpp` still exists; or
* its entry remains in `CMakeLists.txt`.

After changing the implementation layout, explicitly search for remaining references to the removed file name.

For example:

```bash
git grep 'Foo.cpp'
```

or an equivalent repository-wide search.

There should be no build-system reference to the deleted source file.

Also update any other explicit source lists if the repository contains additional build definitions referring to the file.

## Tests

Migrating the tests is part of the migration, not optional cleanup.

Use GoogleTest typed tests to execute the same behavioral contract against a small representative set of template configurations.

Do **not** test the Cartesian product of every possible node and edge-weight type. Choose configurations that exercise meaningful differences.

A typical algorithm test matrix is:

```cpp
using TestTypes = ::testing::Types<
    Config<node, edgeweight>,
    Config<int, float>>;
```

Add further configurations when they expose materially different behavior, for example:

```cpp
Config<int, int>
```

when integral edge weights are supported and numerical behavior needs coverage.

Tests should always include the traditional NetworKit configuration so that the compatibility specialization remains covered.

A useful fixture pattern is:

```cpp
template <typename Config>
class GenericFooGTest : public testing::Test {
protected:
    using GraphT = typename Config::GraphT;
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
};
```

If test configuration requires more than types — for example weighted/unweighted, indexed/unindexed, or directed/undirected graphs — use a small configuration type rather than duplicating complete tests.

Use typed node values explicitly where useful:

```cpp
NodeT{0}
NodeT{1}
static_cast<NodeT>(value)
```

rather than depending on implicit conversion from NetworKit's traditional `node` type.

If readers or other infrastructure have not yet been migrated, reuse an existing generic graph conversion helper instead of duplicating conversion logic in every test.

### GoogleTest assertion style

Use `EXPECT_THAT` extensively for arrays, vectors, ranges, and other container-like results.

Preferred patterns:

| Situation                         | Preferred assertion                                              |
| --------------------------------- | ---------------------------------------------------------------- |
| vector/array equality             | `EXPECT_THAT(actual, testing::ElementsAreArray(expected));`      |
| empty container                   | `EXPECT_THAT(actual, testing::IsEmpty());`                       |
| container size                    | `EXPECT_THAT(actual, testing::SizeIs(expectedSize));`            |
| every element satisfies condition | `EXPECT_THAT(actual, testing::Each(matcher));`                   |
| nested container structure        | nested `Each`, `ElementsAreArray`, or other appropriate matchers |
| scalar equality                   | `EXPECT_EQ(actual, expected);`                                   |
| simple boolean condition          | `EXPECT_TRUE(...)` / `EXPECT_FALSE(...)`                         |
| exception behavior                | `EXPECT_THROW(...)`                                              |
| floating-point scalar             | appropriate floating-point GoogleTest assertion                  |

In particular, do not write:

```cpp
EXPECT_EQ(actualVector, expectedVector);
```

when a suitable gMock matcher provides a clearer container comparison.

Prefer:

```cpp
EXPECT_THAT(actualVector, testing::ElementsAreArray(expectedVector));
```

Likewise, replace manual loops whose sole purpose is asserting the same property for every element:

```cpp
for (...) {
    EXPECT_FALSE(flags[i]);
}
```

with:

```cpp
EXPECT_THAT(flags, testing::Each(false));
```

when this makes the test shorter and clearer.

However, **do not mechanically replace short and expressive `EXPECT_EQ`, `EXPECT_TRUE`, or `EXPECT_FALSE` assertions with `EXPECT_THAT`.**

For example, keep:

```cpp
EXPECT_EQ(test.getDistance(NodeT{0}, NodeT{3}), expectedDistance);
EXPECT_TRUE(test.isPlanar());
EXPECT_FALSE(test.isNodeInNegativeCycle(NodeT{2}));
```

The goal is expressive tests, not maximum matcher usage.

## Test helpers

When many assertions reconstruct a complete result from a scalar API, create a small fixture helper and compare the resulting structure.

For example:

```cpp
std::vector<std::vector<DistanceT>>
distancesFromGetDistance(const AlgorithmT &algorithm, NodeT numberOfNodes) {
    // ...
}
```

followed by:

```cpp
EXPECT_THAT(
    distancesFromGetDistance(test, numberOfNodes),
    testing::ElementsAreArray(expectedDistances));
```

This is preferable to deeply nested loops containing repeated scalar assertions when the actual contract being tested is a complete vector or matrix.

Keep helpers local to the fixture or test translation unit unless multiple test suites genuinely need them.

## Compatibility checks

Before considering the migration complete, verify both interfaces:

```cpp
GenericFoo<SomeGraphType> genericAlgorithm(graph);
Foo legacyAlgorithm(defaultGraph);
```

The legacy class name must continue to compile and behave as before.

Do not require users of the traditional `Graph` API to add template arguments as a consequence of the migration.

If the class participates in inheritance, check inherited APIs with non-default types as well. Generic derived classes can expose signedness, conversion, or override issues in their base classes that the default specialization did not reveal.

Fix such issues at the appropriate abstraction level rather than introducing local workarounds.

## Build and dependency hygiene

Before making the migration, inspect:

* the current `.hpp` and `.cpp` implementation;
* the corresponding `CMakeLists.txt`;
* the current tests;
* already-templated base classes;
* nearby migrated algorithms/data structures;
* direct users of the old class;
* prerequisite migration PRs or commits.

Rebase or update the branch before starting if it depends on another migration. The diff should contain only changes required for the current component.

After moving the implementation into template-visible headers:

```text
1. Delete the original implementation .cpp file.
2. Remove that .cpp file from the corresponding CMakeLists.txt.
3. Search the repository for any remaining references to the deleted file.
4. Build the affected target/module.
5. Run the complete affected GoogleTest suite.
6. Run formatting and repository lint checks.
7. Check for compiler warnings, especially sign comparisons and conversions.
8. Check that the legacy alias still compiles.
9. Inspect the diff for unrelated changes.
10. For performance-sensitive code, compare the default specialization with the pre-migration implementation.
```

Do not declare success merely because the new templated specialization compiles.

A successful build while an obsolete `.cpp` file or stale CMake entry remains does **not** count as a complete migration.

## Performance-sensitive migrations

Templating is partly intended to enable lighter-weight and faster graph representations.

Therefore avoid changes that erase those advantages, such as repeatedly converting `NodeT` back to the traditional 64-bit `node` type inside hot loops.

Prefer containers parameterized by the actual semantic type:

```cpp
std::vector<NodeT>
std::queue<NodeT>
```

rather than:

```cpp
std::vector<node>
std::queue<node>
```

when those objects store node IDs.

For helper data structures, avoid conflating an external value type with an internal container index type.

For hot or fundamental code, benchmark the traditional/default specialization before and after the migration. The migration should normally have no meaningful regression for existing users.

## Scope control

Do not combine the migration with unrelated refactoring merely because templating exposes opportunities for cleanup.

Make small correctness changes when they are necessary for generic instantiation, such as:

* fixing signed/unsigned comparisons;
* replacing a fixed sentinel with a type-aware sentinel;
* initializing previously implicit state required by generic instantiations;
* removing the original `.cpp` implementation;
* removing its `CMakeLists.txt` entry;
* adding missing direct includes.

For larger unrelated improvements, leave the existing behavior intact and handle them separately.

## Final review checklist

Before finishing, confirm all of the following:

* The implementation accepts the intended generic graph/data types.
* Graph-dependent node IDs use `NodeT`.
* Graph-dependent edge weights use `EdgeWeightT` where semantically appropriate.
* Result types that intentionally differ from `EdgeWeightT` are explicit.
* Type-aware sentinels are used.
* The traditional public class name remains available through an alias.
* No unnecessary conversion back to fixed NetworKit types occurs in hot paths.
* Template implementations are visible from headers.
* **The original source `.cpp` implementation has been deleted.**
* **The deleted `.cpp` has also been removed from the corresponding `CMakeLists.txt`.**
* **No stale build-system reference to the deleted `.cpp` remains.**
* Tests are typed over representative configurations.
* The default legacy configuration is included in the test matrix.
* Container and vector comparisons preferentially use `EXPECT_THAT`.
* Short scalar `EXPECT_EQ` and boolean `EXPECT_TRUE`/`EXPECT_FALSE` assertions remain simple.
* The affected tests pass.
* Formatting and compiler warnings are clean.
* The diff contains no accidental unrelated changes.
* Performance-sensitive default specializations show no meaningful regression.

## Expected completion report

When performing a migration using this skill, finish by reporting:

* the new templated class/data-structure name;
* the preserved legacy alias;
* which types were generalized;
* which types intentionally remained fixed and why;
* whether the implementation stayed in the primary header or moved to an `Impl.hpp`;
* the original `.cpp` file that was deleted;
* the corresponding `CMakeLists.txt` entry that was removed;
* the representative typed-test configurations;
* tests and checks executed;
* benchmark results when performance validation was warranted;
* any remaining migration dependency or limitation.

Do not claim the migration is complete if:

* only the implementation was templated while the corresponding tests remain fixed to the legacy `Graph` type;
* the original implementation `.cpp` still exists; or
* the deleted `.cpp` remains referenced by CMake or another build-system source list.
