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

target_link_libraries(networkit_sparsification
    PRIVATE
        OpenMP::OpenMP_CXX
    PUBLIC
        networkit_auxiliary
        networkit_edgescores
        networkit_graph
    )

add_subdirectory(test)

