networkit_add_module(geometric
    HyperbolicSpace.cpp
    )

networkit_module_link_modules(geometric
    auxiliary)

add_subdirectory(test)

