networkit_add_module(correlation
    Assortativity.cpp
    )

networkit_module_link_modules(correlation
    graph structures)

add_subdirectory(test)