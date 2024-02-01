networkit_add_module(bipartite
        Bipartite.cpp
    )

networkit_module_link_modules(graph
        structures
    )

add_subdirectory(test)