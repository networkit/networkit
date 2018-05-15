networkit_add_module(components
    BiconnectedComponents.cpp
    ConnectedComponents.cpp
    DynConnectedComponents.cpp
    DynWeaklyConnectedComponents.cpp
    ParallelConnectedComponents.cpp
    StronglyConnectedComponents.cpp
    WeaklyConnectedComponents.cpp
    )

target_link_libraries(networkit_components
    PUBLIC
        networkit_auxiliary
        networkit_base
        networkit_coarsening
        networkit_distance
        networkit_dynamics
        networkit_graph
        networkit_structures
    )

add_subdirectory(test)

