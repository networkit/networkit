networkit_add_module(embedding
    AliasSampler.cpp
    BiasedRandomWalk.cpp
    LearnEmbedding.cpp
    Node2Vec.cpp
    )

networkit_module_link_modules(embedding
        auxiliary graph)

add_subdirectory(test)

