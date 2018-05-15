networkit_add_module(matching
    LocalMaxMatcher.cpp
    Matcher.cpp
    Matching.cpp
    PathGrowingMatcher.cpp
    )

target_link_libraries(networkit_matching
    PUBLIC
        networkit_auxiliary
        networkit_base
        networkit_graph
        networkit_structures
    )

add_subdirectory(test)

