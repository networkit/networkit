networkit_add_test(graph GraphBuilderAutoCompleteGTest auxiliary)
networkit_add_test(graph GraphBuilderDirectSwapGTest auxiliary)
networkit_add_test(graph GraphGTest
    auxiliary dyn_distance io generators)
networkit_add_test(graph GraphToolsGTest)
networkit_add_test(graph SpanningGTest io)

networkit_add_benchmark(graph Graph2Benchmark)
networkit_add_benchmark(graph GraphBenchmark auxiliary)

