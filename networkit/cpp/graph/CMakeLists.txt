networkit_add_module(graph
    Graph.cpp
    GraphBuilder.cpp
    GraphTools.cpp
    KruskalMSF.cpp
    RandomMaximumSpanningForest.cpp
    RandomSpanningForest.cpp
    Sampling.cpp
    SpanningForest.cpp
    UnionMaximumSpanningForest.cpp
    )

networkit_module_link_libraries(graph
    PRIVATE OpenMP::OpenMP_CXX)

networkit_module_link_modules(graph
    auxiliary base components structures)

add_subdirectory(test)

