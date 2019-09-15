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
#include <tlx/define/likely.hpp>

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
	if (!value) {
		buffer[0] = 1;
		return 1;
	}

	int dataBytes;

	if (TLX_UNLIKELY(value >= (1llu << 56))) {
		// We have more than 56 bits of data, i.e. we need 8 data bytes.
		buffer[0] = 0;
		dataBytes = 8;

	} else {
		const auto bits = 64 - tlx::clz(value);
		dataBytes = (bits - 1) / 7;

		// number of bytes is encode by the number of trailing zeros
		buffer[0] = static_cast<uint8_t>(1 << dataBytes);

		// put the 8 - (dataBytes - 1) least significant bits into the remainder of the buffer byte
		buffer[0] |= static_cast<uint8_t>(value << (dataBytes + 1));
		value >>= 7 - dataBytes;
	}

	for (uint64_t i = 0; i < dataBytes; i++) {
		buffer[i + 1] = static_cast<uint8_t>(value & 0xFF);
		value >>= 8;
	}

	return dataBytes + 1;
}

/// Converts a serialized varint into an uint64_t value and returns the number of bytes consumed.
inline size_t varIntDecode(const uint8_t *data, uint64_t &value) noexcept {
	int n = 8;
	int bitsRecovered = 0;
	uint64_t decoded = 0;

	if (data[0]) {
		n = tlx::ffs(data[0]) - 1;
		decoded = data[0] >> (n + 1);
		bitsRecovered = 7 - n;
	}

	for (int i = 0; i < n; i++) {
		decoded |= static_cast<uint64_t>(data[i + 1]) << bitsRecovered;
		bitsRecovered += 8;
	}

	value = decoded;
	return n + 1;
}

} // namespace nkbg
} // namespace NetworKit

#endif
