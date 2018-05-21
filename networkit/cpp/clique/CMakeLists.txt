networkit_add_module(clique
    MaxClique.cpp
    MaximalCliques.cpp
    )

networkit_module_link_modules(clique
    auxiliary base centrality graph)

add_subdirectory(test)

