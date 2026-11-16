set(NETWORKIT_CENTRALITY_GPU_SOURCES "")
if(NETWORKIT_WITH_CUDA)
    list(APPEND NETWORKIT_CENTRALITY_GPU_SOURCES gpu/CudaPageRank.cpp)
endif()
if(NETWORKIT_WITH_METAL)
    list(APPEND NETWORKIT_CENTRALITY_GPU_SOURCES gpu/MetalPageRank.cpp)
endif()

networkit_add_module(centrality
    ${NETWORKIT_CENTRALITY_GPU_SOURCES}
    ApproxBetweenness.cpp
    ApproxCloseness.cpp
    ApproxElectricalCloseness.cpp
    ApproxGroupBetweenness.cpp
    ApproxSpanningEdge.cpp
    Betweenness.cpp
    Centrality.cpp
    Closeness.cpp
    ComplexPaths.cpp
    CoreDecomposition.cpp
    DegreeCentrality.cpp
    DynApproxBetweenness.cpp
    DynBetweenness.cpp
    DynBetweennessOneNode.cpp
    DynKatzCentrality.cpp
    DynTopHarmonicCloseness.cpp
    EigenvectorCentrality.cpp
    EstimateBetweenness.cpp
    ForestCentrality.cpp
    GedWalk.cpp
    GPUPageRank.cpp
    GroupCloseness.cpp
    GroupClosenessLocalSearch.cpp
    GroupClosenessGrowShrink.cpp
    GroupClosenessGrowShrinkImpl.cpp
    GroupClosenessLocalSwaps.cpp
    GroupDegree.cpp
    GroupHarmonicCloseness.cpp
    HarmonicCloseness.cpp
    KadabraBetweenness.cpp
    KatzCentrality.cpp
    KPathCentrality.cpp
    LaplacianCentrality.cpp
    LocalClusteringCoefficient.cpp
    LocalPartitionCoverage.cpp
    LocalSquareClusteringCoefficient.cpp
    PageRank.cpp
    PermanenceCentrality.cpp
    Sfigality.cpp
    SpanningEdgeCentrality.cpp
    TopCloseness.cpp
    TopHarmonicCloseness.cpp
    )

networkit_module_link_modules(centrality
        algebraic
        auxiliary
        components
        distance
        dyn_distance
        dynamics
        graph
        linkprediction
        reachability
        structures
    )

add_subdirectory(test)

# GPUPageRank.cpp calls only the backends that are built.
set(NETWORKIT_GPU_DEFINITIONS "")
set(NETWORKIT_GPU_LIBRARIES "")

if(NETWORKIT_WITH_CUDA)
    list(APPEND NETWORKIT_GPU_DEFINITIONS NETWORKIT_WITH_CUDA)
    # nvcc compiles the kernels on their own, without the flags meant for the C++ compiler.
    add_library(networkit_centrality_cuda STATIC gpu/CudaPageRankKernels.cu)
    set_target_properties(networkit_centrality_cuda PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        CUDA_RUNTIME_LIBRARY Shared)
    target_link_libraries(networkit_centrality_cuda PUBLIC CUDA::cudart)
    list(APPEND NETWORKIT_GPU_LIBRARIES networkit_centrality_cuda)
endif()

if(NETWORKIT_WITH_METAL)
    list(APPEND NETWORKIT_GPU_DEFINITIONS NETWORKIT_WITH_METAL)

    # The kernels are compiled into a Metal library, whose bytes MetalPageRank.cpp includes, so
    # that there is no file to install and find at runtime.
    set(NETWORKIT_METAL_KERNELS_DIR ${CMAKE_CURRENT_BINARY_DIR}/gpu)
    set(NETWORKIT_METAL_FLAGS "")
    if(CMAKE_OSX_DEPLOYMENT_TARGET)
        list(APPEND NETWORKIT_METAL_FLAGS -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET})
    endif()
    add_custom_command(
        OUTPUT ${NETWORKIT_METAL_KERNELS_DIR}/PageRankMetallib.hpp
        BYPRODUCTS ${NETWORKIT_METAL_KERNELS_DIR}/PageRank.metallib
        COMMAND ${CMAKE_COMMAND} -E make_directory ${NETWORKIT_METAL_KERNELS_DIR}
        COMMAND xcrun -sdk macosx metal ${NETWORKIT_METAL_FLAGS}
                -o ${NETWORKIT_METAL_KERNELS_DIR}/PageRank.metallib
                ${CMAKE_CURRENT_SOURCE_DIR}/gpu/PageRank.metal
        COMMAND ${CMAKE_COMMAND}
                -DINPUT=${NETWORKIT_METAL_KERNELS_DIR}/PageRank.metallib
                -DOUTPUT=${NETWORKIT_METAL_KERNELS_DIR}/PageRankMetallib.hpp
                -DNAME=pageRankMetallib
                -P ${CMAKE_CURRENT_SOURCE_DIR}/gpu/EmbedFile.cmake
        DEPENDS gpu/PageRank.metal gpu/EmbedFile.cmake
        COMMENT "Compiling the Metal kernels of PageRank"
        VERBATIM)
    add_custom_target(networkit_centrality_metal
        DEPENDS ${NETWORKIT_METAL_KERNELS_DIR}/PageRankMetallib.hpp)
    if(NETWORKIT_MONOLITH)
        add_dependencies(networkit networkit_centrality_metal)
    else()
        add_dependencies(networkit_centrality networkit_centrality_metal)
    endif()

    # Only MetalPageRank.cpp sees the metal-cpp headers, which are not warning-free.
    set_source_files_properties(gpu/MetalPageRank.cpp
        DIRECTORY ${PROJECT_SOURCE_DIR}
        PROPERTIES
            COMPILE_OPTIONS "-isystem;${NETWORKIT_METAL_CPP_PATH};-I${NETWORKIT_METAL_KERNELS_DIR}"
            OBJECT_DEPENDS ${NETWORKIT_METAL_KERNELS_DIR}/PageRankMetallib.hpp)
    list(APPEND NETWORKIT_GPU_LIBRARIES "-framework Foundation" "-framework Metal")
endif()

set_source_files_properties(GPUPageRank.cpp
    DIRECTORY ${PROJECT_SOURCE_DIR}
    PROPERTIES COMPILE_DEFINITIONS "${NETWORKIT_GPU_DEFINITIONS}")
if(NETWORKIT_GPU_LIBRARIES)
    if(NETWORKIT_MONOLITH)
        target_link_libraries(networkit PRIVATE ${NETWORKIT_GPU_LIBRARIES})
    else()
        networkit_module_link_libraries(centrality PRIVATE ${NETWORKIT_GPU_LIBRARIES})
    endif()
endif()
