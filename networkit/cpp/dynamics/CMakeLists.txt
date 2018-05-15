networkit_add_module(dynamics
    AffectedNodes.cpp
    DGSStreamParser.cpp
    DGSWriter.cpp
    GraphDifference.cpp
    GraphEvent.cpp
    GraphEventHandler.cpp
    GraphEventProxy.cpp
    GraphUpdater.cpp
    )

target_link_libraries(networkit_dynamics
    PRIVATE
        OpenMP::OpenMP_CXX
    PUBLIC
        networkit_auxiliary
        networkit_base
        networkit_distance
        networkit_graph
    )

add_subdirectory(test)

