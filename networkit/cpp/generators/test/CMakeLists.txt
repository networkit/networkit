networkit_add_test(generators GeneratorsGTest
    auxiliary community dynamics global io viz)

networkit_add_test(generators ErdosRenyiEnumeratorGTest)

networkit_add_benchmark(generators GeneratorsBenchmark
        auxiliary graph)
