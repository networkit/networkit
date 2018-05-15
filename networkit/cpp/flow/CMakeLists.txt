networkit_add_module(flow
    EdmondsKarp.cpp
    )

target_link_libraries(networkit_flow
    PUBLIC
        networkit_graph
    )

add_subdirectory(test)

