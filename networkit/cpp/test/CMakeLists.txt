networkit_add_module(test
        BasicsBenchmarkGTest.cpp
        ParametrizedGTest.cpp
        )

target_link_libraries(networkit_test
    PRIVATE
        GTest::GTest
        GTest::Main
    PUBLIC
        networkit_auxiliary
    )
