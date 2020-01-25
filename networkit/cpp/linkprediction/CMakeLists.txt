networkit_add_module(linkprediction
    AdamicAdarIndex.cpp
    AdjustedRandIndex.cpp
    AlgebraicDistanceIndex.cpp
    EvaluationMetric.cpp
    KatzIndex.cpp
    LinkPredictor.cpp
    LinkThresholder.cpp
    MissingLinksFinder.cpp
    NeighborhoodUtility.cpp
    PrecisionRecallMetric.cpp
    PredictionsSorter.cpp
    ROCMetric.cpp
    RandomLinkSampler.cpp
    ResourceAllocationIndex.cpp
    SameCommunityIndex.cpp
    )

networkit_module_link_modules(linkprediction
    auxiliary community graph structures)

add_subdirectory(test)

