networkit_add_module(geometric
    HyperbolicSpace.cpp
    )

target_link_libraries(networkit_geometric
    PUBLIC
        networkit_auxiliary
    )

add_subdirectory(test)

