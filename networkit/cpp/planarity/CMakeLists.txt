networkit_add_module(planarity
    LeftRightPlanarityTest.cpp
    )

networkit_module_link_modules(planarity
    auxiliary coarsening dynamics graph structures)

add_subdirectory(test)

