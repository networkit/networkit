#ifndef NETWORKIT_BINARY_GRAPH_HPP
#define NETWORKIT_BINARY_GRAPH_HPP

namespace nkbg {
struct Header {
	char magic[8];
	uint64_t checksum;
	uint64_t features;
	uint64_t nodes;
	uint64_t chunks;
	uint64_t offsetBaseData;
	uint64_t offsetAdjLists;
	uint64_t offsetAdjTranspose;
	uint64_t offsetWeightLists;
	uint64_t offsetWeightTranspose;
};

enum WEIGHT_FORMAT {
	NONE = 0,
	VARINT = 1,
	SIGNED_VARINT = 2,
	DOUBLE = 3,
	FLOAT = 4
};

static constexpr uint8_t DELETED_BIT = 0x1; // bit 0
static constexpr uint64_t DIR_MASK = 0x1; // bit 0
static constexpr uint64_t WGHT_MASK = 0xE; //bit 1-3
static constexpr uint64_t WGHT_SHIFT = 0x1;

}
#endif
