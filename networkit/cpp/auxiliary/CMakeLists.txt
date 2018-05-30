networkit_add_module(auxiliary
    BloomFilter.cpp
    BucketPQ.cpp
    Log.cpp
    NumericTools.cpp
    Parallelism.cpp
    PrioQueueForInts.cpp
    Random.cpp
    SignalHandling.cpp
    Timer.cpp
    )

networkit_module_link_libraries(auxiliary
    PRIVATE
        OpenMP::OpenMP_CXX
    )

add_subdirectory(test)

