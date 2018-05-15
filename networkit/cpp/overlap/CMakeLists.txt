networkit_add_module(overlap
    HashingOverlapper.cpp
    Overlapper.cpp
    )

target_link_libraries(networkit_overlap
    PUBLIC
        networkit_auxiliary
        networkit_graph
        networkit_structures
    )

add_subdirectory(test)

