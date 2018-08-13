networkit_add_module(io
    BinaryEdgeListPartitionReader.cpp
    BinaryEdgeListPartitionWriter.cpp
    BinaryPartitionReader.cpp
    BinaryPartitionWriter.cpp
    CoverReader.cpp
    CoverWriter.cpp
    DGSReader.cpp
    DibapGraphReader.cpp
    DotGraphWriter.cpp
    DotPartitionWriter.cpp
    DynamicGraphReader.cpp
    EdgeListCoverReader.cpp
    EdgeListPartitionReader.cpp
    EdgeListReader.cpp
    EdgeListWriter.cpp
    GMLGraphReader.cpp
    GMLGraphWriter.cpp
    GraphIO.cpp
    GraphToolBinaryReader.cpp
    GraphToolBinaryWriter.cpp
    KONECTGraphReader.cpp
    LineFileReader.cpp
    MemoryMappedFile.cpp
    METISGraphReader.cpp
    METISGraphWriter.cpp
    METISParser.cpp
    MatrixMarketReader.cpp
    PartitionReader.cpp
    PartitionWriter.cpp
    RasterReader.cpp
    SNAPEdgeListPartitionReader.cpp
    SNAPGraphReader.cpp
    SNAPGraphWriter.cpp
    ThrillGraphBinaryReader.cpp
    ThrillGraphBinaryWriter.cpp
    )

networkit_module_link_modules(io
    algebraic auxiliary dynamics graph structures)

add_subdirectory(test)

