networkit_add_module(clique
    MaximalCliques.cpp
    )

networkit_module_link_modules(clique
    auxiliary base centrality graph)

add_subdirectory(test)

