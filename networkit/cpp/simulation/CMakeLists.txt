networkit_add_module(simulation
    EpidemicSimulationSEIR.cpp
    )

target_link_libraries(networkit_simulation
    PUBLIC
        networkit_auxiliary
        networkit_base
        networkit_graph
    )

