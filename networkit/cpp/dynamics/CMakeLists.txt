networkit_add_module(dynamics
    DGSStreamParser.cpp
    DGSWriter.cpp
    GraphDifference.cpp
    GraphEvent.cpp
    GraphEventProxy.cpp
    GraphUpdater.cpp
    )

networkit_module_link_modules(dynamics
        auxiliary graph)

add_subdirectory(test)

