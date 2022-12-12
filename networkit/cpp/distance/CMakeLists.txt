networkit_add_module(distance
    AffectedNodes.cpp
    APSP.cpp
    AdamicAdarDistance.cpp
    AlgebraicDistance.cpp
    BFS.cpp
    BidirectionalBFS.cpp
    BidirectionalDijkstra.cpp
    CommuteTimeDistance.cpp
    Diameter.cpp
    Dijkstra.cpp
    Eccentricity.cpp
    EffectiveDiameter.cpp
    EffectiveDiameterApproximation.cpp
    GraphDistance.cpp
    HopPlotApproximation.cpp
    IncompleteDijkstra.cpp
    JaccardDistance.cpp
    MultiTargetBFS.cpp
    MultiTargetDijkstra.cpp
    NeighborhoodFunction.cpp
    NeighborhoodFunctionApproximation.cpp
    NeighborhoodFunctionHeuristic.cpp
    PrunedLandmarkLabeling.cpp
    ReverseBFS.cpp
    SPSP.cpp
    SSSP.cpp
    Volume.cpp
    )

networkit_module_link_modules(distance
        algebraic auxiliary components graph structures)

networkit_add_module(dyn_distance
    DynAPSP.cpp
    DynBFS.cpp
    DynDijkstra.cpp
    DynSSSP.cpp
    )

networkit_module_link_modules(dyn_distance
        distance dynamics)

add_subdirectory(test)
