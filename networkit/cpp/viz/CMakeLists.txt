networkit_add_module(viz
    FruchtermanReingold.cpp
    Layouter.cpp
    MaxentStress.cpp
    MultilevelLayouter.cpp
    PivotMDS.cpp
    PostscriptWriter.cpp
    )

networkit_module_link_modules(viz
    algebraic auxiliary coarsening community components distance graph io structures)

add_subdirectory(test)
