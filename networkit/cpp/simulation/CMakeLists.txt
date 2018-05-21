networkit_add_module(simulation
    EpidemicSimulationSEIR.cpp
    )

networkit_module_link_modules(simulation
    auxiliary base graph)

