networkit_add_module(clique
    MaxClique.cpp
    MaximalCliques.cpp
    )

target_link_libraries(networkit_clique
    PUBLIC
        networkit_auxiliary
        networkit_base
        networkit_centrality
        networkit_graph
    )

add_subdirectory(test)

