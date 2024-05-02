networkit_add_module(numerics
    MultiLevelSetup.cpp
    LevelHierarchy.cpp
    SolverLamg.cpp
    Lamg.cpp
)

networkit_module_link_modules(numerics
    algebraic)

add_subdirectory(test)

