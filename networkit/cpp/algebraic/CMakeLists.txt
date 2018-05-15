networkit_add_module(algebraic
    CSRMatrix.cpp
    DenseMatrix.cpp
    DynamicMatrix.cpp
    Vector.cpp
    )

target_link_libraries(networkit_algebraic
    PRIVATE
        OpenMP::OpenMP_CXX
    PUBLIC
        networkit_auxiliary
        networkit_graph
    )

add_subdirectory(algorithms)
add_subdirectory(test)

