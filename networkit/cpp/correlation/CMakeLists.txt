networkit_add_module(correlation
    Assortativity.cpp
    )

target_link_libraries(networkit_correlation
    PUBLIC
        networkit_base
        networkit_graph
        networkit_structures
    )

