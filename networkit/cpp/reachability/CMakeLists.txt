networkit_add_module(reachability
    AllSimplePaths.cpp
    ReachableNodes.cpp
    )

networkit_module_link_modules(reachability
    components
    graph
    )

add_subdirectory(test)
