networkit_add_module(structures
    Cover.cpp
    Partition.cpp
    UnionFind.cpp
    )

networkit_module_link_libraries(structures
    PRIVATE OpenMP::OpenMP_CXX)


add_subdirectory(test)

