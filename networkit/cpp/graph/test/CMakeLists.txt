networkit_add_test(graph GraphBuilderAutoCompleteGTest auxiliary)
networkit_add_test(graph GraphGTest
    auxiliary dyn_distance io generators)
networkit_add_test(graph GraphToolsGTest generators io)
networkit_add_test(graph TraversalGTest generators)
networkit_add_test(graph SpanningGTest io)
networkit_add_test(graph TopologicalSortGTest)
networkit_add_test(graph AttributeGTest graph)
networkit_add_test(graph HypergraphGTest generators graph)
networkit_add_test(graph HypergraphToolsGTest generators graph)

networkit_add_benchmark(graph Graph2Benchmark)
networkit_add_benchmark(graph GraphBenchmark auxiliary)
