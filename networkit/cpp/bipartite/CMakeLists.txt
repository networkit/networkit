networkit_add_module(bipartite BFSBipartiteCheck.cpp)

networkit_module_link_modules(bipartite graph)

add_subdirectory(test)

