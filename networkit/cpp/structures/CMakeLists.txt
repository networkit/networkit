networkit_add_module(structures
    Cover.cpp
    Partition.cpp
    UnionFind.cpp
    )

target_link_libraries(networkit_structures
    PRIVATE
        OpenMP::OpenMP_CXX
    PUBLIC
        networkit_graph
    )

add_subdirectory(test)

