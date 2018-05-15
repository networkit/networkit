networkit_add_module(independentset
    IndependentSetFinder.cpp
    Luby.cpp
    )

target_link_libraries(networkit_independentset
    PUBLIC
        networkit_auxiliary
        networkit_graph
    )

add_subdirectory(test)

