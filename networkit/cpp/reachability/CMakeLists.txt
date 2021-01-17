networkit_add_module(reachability
    ReachableNodes.cpp
    )

networkit_module_link_modules(reachability
    base
    components
    graph
    )

add_subdirectory(test)
