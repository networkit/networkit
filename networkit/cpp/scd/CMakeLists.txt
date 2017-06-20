networkit_add_module(scd
    ApproximatePageRank.cpp
    CliqueDetect.cpp
    CombinedSCD.cpp
    GCE.cpp
    LFMLocal.cpp
    PageRankNibble.cpp
    SCDGroundTruthComparison.cpp
    SelectiveCommunityDetector.cpp
    SetConductance.cpp
    )

networkit_module_link_modules(scd
    auxiliary clique graph)

add_subdirectory(test)
