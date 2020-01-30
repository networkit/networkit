networkit_add_module(overlap
    HashingOverlapper.cpp
    )

networkit_module_link_modules(overlap
    auxiliary graph structures)

add_subdirectory(test)

