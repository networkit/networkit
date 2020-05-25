networkit_add_module(distance
    AffectedNodes.cpp
    APSP.cpp
    AdamicAdarDistance.cpp
    AlgebraicDistance.cpp
    AllSimplePaths.cpp
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
    NeighborhoodFunction.cpp
    NeighborhoodFunctionApproximation.cpp
    NeighborhoodFunctionHeuristic.cpp
    ReverseBFS.cpp
    SSSP.cpp
    Volume.cpp
    )

networkit_module_link_modules(distance
        algebraic auxiliary base components graph structures)

networkit_add_module(dyn_distance
    DynAPSP.cpp
    DynBFS.cpp
    DynDijkstra.cpp
    DynSSSP.cpp
    )

networkit_module_link_modules(dyn_distance
        distance dynamics)

add_subdirectory(test)
