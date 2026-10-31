networkit_add_test(io IOGTest
        auxiliary community generators structures)

networkit_add_test(io IOTestUtilsGTest)

networkit_add_test(io METISGraphIOGTest
        generators)

networkit_add_test(io METISParserGTest)

networkit_add_test(io MemoryMappedFileGTest)

networkit_add_benchmark(io IOBenchmark geometric)

