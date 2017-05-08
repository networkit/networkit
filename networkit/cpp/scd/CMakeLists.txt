networkit_add_module(scd
    ApproximatePageRank.cpp
    GCE.cpp
    PageRankNibble.cpp
    SCDGroundTruthComparison.cpp
    )

networkit_module_link_modules(scd
    auxiliary graph)

add_subdirectory(test)
