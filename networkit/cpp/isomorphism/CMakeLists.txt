networkit_add_module(isomorphism
    SearchGraph.cpp
    SubgraphIsomorphism.cpp
    )

networkit_module_link_modules(isomorphism
    auxiliary graph)

add_subdirectory(test)