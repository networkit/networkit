networkit_add_module(global
    ClusteringCoefficient.cpp
    GlobalClusteringCoefficient.cpp
    )

networkit_module_link_modules(global
    auxiliary centrality graph)

add_subdirectory(test)

