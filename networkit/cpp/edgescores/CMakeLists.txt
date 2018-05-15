networkit_add_module(edgescores
    ChibaNishizekiQuadrangleEdgeScore.cpp
    ChibaNishizekiTriangleEdgeScore.cpp
    EdgeScore.cpp
    EdgeScoreAsWeight.cpp
    EdgeScoreBlender.cpp
    EdgeScoreLinearizer.cpp
    EdgeScoreNormalizer.cpp
    GeometricMeanScore.cpp
    PrefixJaccardScore.cpp
    TriangleEdgeScore.cpp
    )

target_link_libraries(networkit_edgescores
    PRIVATE
        OpenMP::OpenMP_CXX
    PUBLIC
        networkit_auxiliary
        networkit_base
        networkit_graph
    )

add_subdirectory(test)

