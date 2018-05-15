networkit_add_module(coarsening
    ClusteringProjector.cpp
    GraphCoarsening.cpp
    MatchingCoarsening.cpp
    ParallelPartitionCoarsening.cpp
    )

target_link_libraries(networkit_coarsening
    PRIVATE
        OpenMP::OpenMP_CXX
    PUBLIC
        networkit_auxiliary
        networkit_base
        networkit_graph
        networkit_matching
        networkit_structures
    )

add_subdirectory(test)

