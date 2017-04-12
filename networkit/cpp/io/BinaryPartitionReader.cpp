#include "BinaryPartitionReader.h"
#include <fstream>


NetworKit::BinaryPartitionReader::BinaryPartitionReader(uint8_t width) : width(width) {
	if (width != 4 && width != 8) {
		throw std::runtime_error("Only 4 and 8 are supported widths");
	}
}

namespace {
	template <uint8_t width>
	NetworKit::Partition readPartition(const std::string& path) {
		using read_t = typename std::conditional<width == 4, uint32_t, uint64_t>::type;
		static_assert(sizeof(read_t) == width, "Error, width is not the width of read_node_t");

		std::ifstream is(path, std::ios_base::in | std::ios_base::binary);

		if (!is) {
			throw std::runtime_error("Error: partition file couldn't be opened");
		}

		is.seekg(0, std::ios_base::end);
		NetworKit::count length = is.tellg() / width;
		is.seekg(0);

		NetworKit::Partition zeta(length);
		
		read_t p = 0;
		for (read_t u = 0; u < length; ++u) {
			is.read(reinterpret_cast<char*>(&p), width);

			if (p >= zeta.upperBound()) {
				zeta.setUpperBound(p + 1);
			}

			zeta[u] = p;
		}

		return zeta;
	};
};

NetworKit::Partition NetworKit::BinaryPartitionReader::read(const std::string& path) {
	switch (width) {
		case 4:
			return readPartition<4>(path);
		default:
			return readPartition<8>(path);
	};
}
