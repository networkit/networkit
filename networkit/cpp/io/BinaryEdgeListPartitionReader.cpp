#include "BinaryEdgeListPartitionReader.h"
#include <fstream>

NetworKit::BinaryEdgeListPartitionReader::BinaryEdgeListPartitionReader(node firstNode, uint8_t width) : firstNode(firstNode), width(width) {
	if (width != 4 && width != 8) {
		throw std::runtime_error("Error: width must be 4 or 8");
	}
}

namespace {

	template <uint8_t width>
	NetworKit::Partition readPartition(const std::string& path, NetworKit::node firstNode) {
		using read_t = typename std::conditional<width == 4, uint32_t, uint64_t>::type;
		static_assert(sizeof(read_t) == width, "Error, width is not the width of read_node_t");

		std::ifstream is(path, std::ios_base::in | std::ios_base::binary);

		if (!is) {
			throw std::runtime_error("Error: partition file couldn't be opened");
		}

		is.seekg(0, std::ios_base::end);
		NetworKit::count length = is.tellg() / (2 * width);
		is.seekg(0);

		NetworKit::Partition zeta(length);
		
		read_t u = 0, p = 0;

		while (is.good()) {
			is.read(reinterpret_cast<char*>(&u), width);
			is.read(reinterpret_cast<char*>(&p), width);

			if (u < firstNode) {
				throw std::runtime_error("Error: node smaller than the given firstNode found!");
			}

			u -= firstNode;

			if (p >= zeta.upperBound()) {
				zeta.setUpperBound(p + 1);
			}

			while (u >= zeta.numberOfElements()) {
				zeta.extend();
			}

			zeta[u] = p;
		}

		return zeta;
	};
};

NetworKit::Partition NetworKit::BinaryEdgeListPartitionReader::read(const std::string& path) {
	switch (width) {
		case 4:
			return readPartition<4>(path, firstNode);
		default:
			return readPartition<8>(path, firstNode);
	};
}
