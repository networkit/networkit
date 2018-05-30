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

networkit_module_link_libraries(dynamics
    PRIVATE OpenMP::OpenMP_CXX)

networkit_module_link_modules(dynamics
        auxiliary base distance graph)

add_subdirectory(test)

