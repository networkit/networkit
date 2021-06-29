networkit_add_module(scd
    ApproximatePageRank.cpp
    CliqueDetect.cpp
    CombinedSCD.cpp
    GCE.cpp
    LFMLocal.cpp
    LocalT.cpp
    LocalTightnessExpansion.cpp
    PageRankNibble.cpp
    RandomBFS.cpp
    SCDGroundTruthComparison.cpp
    SelectiveCommunityDetector.cpp
    SetConductance.cpp
    TwoPhaseL.cpp
    TCE.cpp
    )

networkit_module_link_modules(scd
    auxiliary clique graph)

add_subdirectory(test)
