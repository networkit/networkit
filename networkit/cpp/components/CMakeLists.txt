networkit_add_module(components
    BiconnectedComponents.cpp
    ConnectedComponents.cpp
    DynConnectedComponents.cpp
    DynWeaklyConnectedComponents.cpp
    ParallelConnectedComponents.cpp
    RandomSpanningForest.cpp
    StronglyConnectedComponents.cpp
    WeaklyConnectedComponents.cpp
    )

networkit_module_link_modules(components
    auxiliary base coarsening dynamics graph structures)

add_subdirectory(test)

