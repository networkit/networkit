networkit_add_module(base
    Algorithm.cpp
    )

target_link_libraries(networkit_base
    PUBLIC
        networkit_auxiliary
        networkit_dynamics
    )

