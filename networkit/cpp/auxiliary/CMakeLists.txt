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

target_link_libraries(networkit_auxiliary
    PRIVATE
        OpenMP::OpenMP_CXX
    )

add_subdirectory(test)

