#ifndef NETWORKIT_BINARY_GRAPH_
#define NETWORKIT_BINARY_GRAPH_

struct Header {
	char magic[8];
	uint64_t checksum;
	uint64_t features; 	
	uint64_t nodes;
	uint64_t chunks;
	uint64_t offsetBaseData;
	uint64_t offsetAdjLists;
	uint64_t offsetAdjTranspose;
	uint64_t offsetWeights;
};

enum WEIGHT_FORMAT {
	VARINT = 0,
	SIGNED_VARINT = 1,
	DOUBLE = 2,
	FLOAT = 3
};

enum MASKS {
	DIR_MASK = 0x1, // bit 0
	WGHT_MASK = 0x6, // bit 1-2
	WGHT_SHIFT = 0x1
};

constexpr auto DELETED_BIT = uint64_t(1) << 52;
constexpr auto SIZE_MASK = (uint64_t(1) << 52) - 1;

#endif
