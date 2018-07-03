networkit_add_module(sparsification
    ChanceCorrectedTriangleScore.cpp
    ForestFireScore.cpp
    GlobalThresholdFilter.cpp
    LocalDegreeScore.cpp
    LocalSimilarityScore.cpp
    MultiscaleScore.cpp
    RandomEdgeScore.cpp
    RandomNodeEdgeScore.cpp
    SCANStructuralSimilarityScore.cpp
    SimmelianOverlapScore.cpp
    SimmelianScore.cpp
    Sparsifiers.cpp
    )

networkit_module_link_modules(sparsification
    auxiliary edgescores graph)

add_subdirectory(test)

