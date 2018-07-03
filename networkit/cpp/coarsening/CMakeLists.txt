networkit_add_module(coarsening
    ClusteringProjector.cpp
    GraphCoarsening.cpp
    MatchingCoarsening.cpp
    ParallelPartitionCoarsening.cpp
    )

networkit_module_link_modules(coarsening
    auxiliary base graph matching structures)

add_subdirectory(test)

