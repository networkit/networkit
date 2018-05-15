networkit_add_module(global
    ClusteringCoefficient.cpp
    GlobalClusteringCoefficient.cpp
    )

target_link_libraries(networkit_global
    PUBLIC
        networkit_auxiliary
        networkit_centrality
        networkit_graph
    )

add_subdirectory(test)

