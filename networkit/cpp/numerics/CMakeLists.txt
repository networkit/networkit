networkit_add_module(numerics
    MultiLevelSetup.cpp
    LevelHierarchy.cpp
)

networkit_module_link_modules(numerics
    algebraic)

add_subdirectory(test)

