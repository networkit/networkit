/*
 * NetworkitBinaryGraph.hpp
 *
 *      Author: Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#ifndef NETWORKIT_BINARY_GRAPH_HPP
#define NETWORKIT_BINARY_GRAPH_HPP

#include <cstdint>
#include <cstddef>

#include <tlx/math/clz.hpp>
#include <tlx/math/ffs.hpp>

namespace NetworKit {
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

/// Serializes value into a buffer and returns the number of bytes written.
inline size_t varIntEncode(uint64_t value, uint8_t *buffer) noexcept {
	uint64_t n;
	if(!value) {
		buffer[0] =1;
		return 1;
	}

	if (value > (uint64_t(1) << 63)) {
		n = 8;
		buffer[0] = 0;
	} else {
		n = (64 - tlx::clz(value)) / 7;
		buffer[0] = (1 << n) | (value << (n + 1));
		value >>= 8 - (n + 1);
	}

	for(uint64_t i = 0; i < n; i++) {
		buffer[i+1] = value & 0xFF;
		value >>= 8;
	}
	return n + 1;
}

/// Converts a serialized varint into an uint64_t value and returns the number of bytes consumed.
inline size_t varIntDecode(const uint8_t *data, uint64_t &value) noexcept {
	int n;
	if(!data[0]) {
		n = 8;
		value = 0;
	} else {
		n = tlx::ffs(data[0]) -1;
		value = data[0] >> (n+1);
	}

	for(int i = 0; i < n; i++) {
		value |= data[i+1] << (8 - (n + 1) + i * 8);
	}
	return n+1;
}

} // namespace nkbg
} // namespace NetworKit

#endif
