networkit_add_module(overlap
    HashingOverlapper.cpp
    Overlapper.cpp
    )

networkit_module_link_libraries(overlap
    auxiliary graph structures)

add_subdirectory(test)

