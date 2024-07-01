networkit_add_module(generators
    BarabasiAlbertGenerator.cpp
    ChungLuGenerator.cpp
    ChungLuGeneratorAlamEtAl.cpp
    ClusteredRandomGraphGenerator.cpp
    ConfigurationModel.cpp
    DorogovtsevMendesGenerator.cpp
    DynamicBarabasiAlbertGenerator.cpp
    DynamicDGSParser.cpp
    DynamicDorogovtsevMendesGenerator.cpp
    DynamicForestFireGenerator.cpp
    DynamicGraphSource.cpp
    DynamicHyperbolicGenerator.cpp
    DynamicPathGenerator.cpp
    DynamicPubWebGenerator.cpp
    EdgeSwitchingMarkovChainGenerator.cpp
    ErdosRenyiGenerator.cpp
    HavelHakimiGenerator.cpp
    HyperbolicGenerator.cpp
    LFRGenerator.cpp
    MocnikGenerator.cpp
    MocnikGeneratorBasic.cpp
    PowerlawDegreeSequence.cpp
    PubWebGenerator.cpp
    RegularRingLatticeGenerator.cpp
    RmatGenerator.cpp
    StaticDegreeSequenceGenerator.cpp
    StochasticBlockmodel.cpp
    WattsStrogatzGenerator.cpp
    )

networkit_module_link_modules(generators
        auxiliary dynamics geometric graph randomization structures)

add_subdirectory(quadtree)
add_subdirectory(test)
