networkit_add_module(randomization
    Curveball.cpp
    CurveballGlobalTradeGenerator.cpp
    CurveballImpl.cpp
    CurveballUniformTradeGenerator.cpp
    GlobalCurveball.cpp
    )

networkit_module_link_modules(randomization
    auxiliary graph)

add_subdirectory(test)

