networkit_add_module(viz
    FruchtermanReingold.cpp
    Layouter.cpp
    MaxentStress.cpp
    MultilevelLayouter.cpp
    PivotMDS.cpp
    PostscriptWriter.cpp
    )

target_link_libraries(networkit_viz
    PRIVATE
        OpenMP::OpenMP_CXX
    PUBLIC
        networkit_algebraic
        networkit_auxiliary
        networkit_coarsening
        networkit_community
        networkit_components
        networkit_distance
        networkit_graph
        networkit_io
        networkit_structures
    )

add_subdirectory(test)

