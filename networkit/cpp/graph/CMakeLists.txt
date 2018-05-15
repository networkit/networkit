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

target_link_libraries(networkit_graph
    PRIVATE
        OpenMP::OpenMP_CXX
    PUBLIC
        networkit_auxiliary
        networkit_base
        networkit_components
        networkit_structures
    )

add_subdirectory(test)

