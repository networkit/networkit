networkit_add_module(isomorphism
    SubgraphIsomorphism.cpp
    )

networkit_module_link_modules(isomorphism
    auxiliary graph)

add_subdirectory(test)