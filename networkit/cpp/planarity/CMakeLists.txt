networkit_add_module(planarity LeftRightPlanarityCheck.cpp)

networkit_module_link_modules(planarity graph)

add_subdirectory(test)

