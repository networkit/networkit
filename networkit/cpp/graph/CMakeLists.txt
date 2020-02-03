networkit_add_module(graph
    Graph.cpp
    GraphBuilder.cpp
    GraphTools.cpp
    KruskalMSF.cpp
    RandomMaximumSpanningForest.cpp
    SpanningForest.cpp
    UnionMaximumSpanningForest.cpp
    )

networkit_module_link_modules(graph
    auxiliary base structures)

add_subdirectory(test)

