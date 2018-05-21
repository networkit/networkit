networkit_add_test(sparsification LocalDegreeGTest)
networkit_add_test(sparsification LocalSimilarityGTest
    edgescores)
networkit_add_test(sparsification MultiscaleBackboneGTest)
networkit_add_test(sparsification SimmelianBackboneGTest
    edgescores)

networkit_add_benchmark(sparsification SparsificationBenchmark
    auxiliary edgescores io)

