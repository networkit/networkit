networkit_add_module(graph
    Graph.cpp
    GraphBuilder.cpp
    GraphTools.cpp
    Hypergraph.cpp
    HypergraphTools.cpp
    KruskalMSF.cpp
    RandomMaximumSpanningForest.cpp
    SpanningForest.cpp
    TopologicalSort.cpp
    UnionMaximumSpanningForest.cpp
)

networkit_module_link_modules(graph
    auxiliary structures)

add_subdirectory(test)
