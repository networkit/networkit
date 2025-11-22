networkit_add_module(planarity LeftRightPlanarityCheck.cpp NewLeftRightPlanarityCheck.cpp)

networkit_module_link_modules(planarity graph)

add_subdirectory(test)

