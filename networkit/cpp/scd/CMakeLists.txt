networkit_add_module(scd
    ApproximatePageRank.cpp
    GCE.cpp
    PageRankNibble.cpp
    SelectiveCommunityDetector.cpp
    )

target_link_libraries(networkit_scd
    PUBLIC
        networkit_auxiliary
        networkit_graph
    )

add_subdirectory(test)

