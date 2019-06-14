NetworkitBinaryGraph format
==========================
This is a documentation of the NetworkitBinaryGraph file format. 
The file starts with a header followed by the base data, adjacency lists and weights respectively. 
The rest of this documentation describes the aforementioned blocks in detail.

Header
---------
The file begins with a header that is structured as follows:
```
struct Header {

    char magic[8];
    uint64_t checksum;
    uint64_t features; 
    uint64_t nodes;
    uint64_t chunks;
    uint64_t offsetBaseData;
    uint64_t offsetAdjLists;
    uint64_t offsetTranspose;
	uint64_t offsetWeights;
};
```
- magic: A constant value used to identify the file format version.
    - The current version is '*nkbg000*' which supports unweighted, undirected graphs.
- checksum: Currently not used
- features: Contains the graph information bitwise
    - Bit 0 : directed or undirected
    - Bit 1-2 : the weight format of the graph. Below are the possible formats:
        - 0 = Weights are unsigned integers
        - 1 = Weights are signed integers
        - 2 = Weights are doubles
        - 3 = Weights are floats
- nodes : The number of nodes the graph has
- chunks: The number of chunks the nodes have been divided in
- offsetBaseData: Offset of base data in the file 
- offsetAdjLists: Offset of the adjacency lists in the file
- offsetTranspose: Currently unused
- offsetWeights: Offset of the weights in the file

All offsets are relative to the beginning of the file.

Base data
------------
```
uint64_t sizeSum[nodes]: The sum of sizes of adjacency vertices with indices <= i
uint64_t firstVertex[chunks-1]: The index of the first vertex of each chunk excluding the first chunk
```
Adjacency lists
-----------------
```
uint64_t offset[chunks-1]: Offset of the file where the adjacency list of 
the firstVertex of each chunk relative to data starts
varint data [...]: Varint encoded adjaceny lists  
```
Weights
--------------------
Weights are currently not supported. 

