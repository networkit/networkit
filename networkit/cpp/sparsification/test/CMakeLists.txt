networkit_add_test(sparsification LocalDegreeGTest)
networkit_add_test(sparsification GlobalThresholdFilterGTest
    generators)
networkit_add_test(sparsification LocalSimilarityGTest
    edgescores)
networkit_add_test(sparsification MultiscaleBackboneGTest)
networkit_add_test(sparsification ForestFireGTest)
networkit_add_test(sparsification RandomNodeEdgeGTest)
networkit_add_test(sparsification TriangleScoreGTest)
networkit_add_test(sparsification SimmelianBackboneGTest
    edgescores)

networkit_add_benchmark(sparsification SparsificationBenchmark
    auxiliary edgescores io)

