networkit_add_module(linkprediction
    AdamicAdarIndex.cpp
    AdjustedRandIndex.cpp
    AlgebraicDistanceIndex.cpp
    CommonNeighborsIndex.cpp
    EvaluationMetric.cpp
    JaccardIndex.cpp
    KatzIndex.cpp
    LinkPredictor.cpp
    LinkThresholder.cpp
    MissingLinksFinder.cpp
    NeighborhoodDistanceIndex.cpp
    NeighborhoodUtility.cpp
    NeighborsMeasureIndex.cpp
    PrecisionRecallMetric.cpp
    PredictionsSorter.cpp
    PreferentialAttachmentIndex.cpp
    ROCMetric.cpp
    RandomLinkSampler.cpp
    ResourceAllocationIndex.cpp
    SameCommunityIndex.cpp
    TotalNeighborsIndex.cpp
    UDegreeIndex.cpp
    VDegreeIndex.cpp
    )

networkit_module_link_modules(linkprediction
    auxiliary community graph structures)

add_subdirectory(test)

