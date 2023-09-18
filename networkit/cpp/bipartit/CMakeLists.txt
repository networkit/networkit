networkit_add_module(bipartit
        Bipartit.cpp
    )

networkit_module_link_modules(graph
        structures
    )

add_subdirectory(test)