networkit_add_module(viz
    MaxentStress.cpp
    PivotMDS.cpp
    PostscriptWriter.cpp
    )

networkit_module_link_modules(viz
    algebraic auxiliary coarsening community components distance graph io structures)

add_subdirectory(test)
