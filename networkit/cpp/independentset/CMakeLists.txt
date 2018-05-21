networkit_add_module(independentset
    IndependentSetFinder.cpp
    Luby.cpp
    )

networkit_module_link_modules(independentset
    auxiliary graph)

add_subdirectory(test)

