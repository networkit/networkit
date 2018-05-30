networkit_add_module(algebraic
    CSRMatrix.cpp
    DenseMatrix.cpp
    DynamicMatrix.cpp
    Vector.cpp
    )

networkit_module_link_libraries(algebraic
    PRIVATE OpenMP::OpenMP_CXX)

networkit_module_link_modules(algebraic
    auxiliary graph)

add_subdirectory(algorithms)
add_subdirectory(test)

