#include "BinaryEdgeListPartitionReader.h"
#include <fstream>

NetworKit::BinaryEdgeListPartitionReader::BinaryEdgeListPartitionReader(node firstNode, uint8_t width) : firstNode(firstNode), width(width) {
	if (width != 4 && width != 8) {
		throw std::runtime_error("Error: width must be 4 or 8");
	}
}

namespace {

	template <uint8_t width>
	NetworKit::Partition readPartition(const std::vector<std::string>& paths, NetworKit::node firstNode) {
		using read_t = typename std::conditional<width == 4, uint32_t, uint64_t>::type;
		static_assert(sizeof(read_t) == width, "Error, width is not the width of read_node_t");

		NetworKit::Partition zeta(0);

		if (!paths.empty()) {
			std::ifstream is;

			NetworKit::count _file_index = 0;

			auto next_input = [&]() {
				is.close();
				if (_file_index < paths.size()) {
					is.open(paths[_file_index++], std::ios_base::in | std::ios_base::binary);
					if (!is) {
						throw std::runtime_error("Error: partition file couldn't be opened");
					}
				}
			};

			next_input();


			read_t u = 0, p = 0;

			while (is.good() && is.is_open()) {
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

				if (is.is_open() && (is.peek() == std::char_traits<char>::eof() || !is.good())) {
					next_input();
				}
			}
		}

		return zeta;
	};
};

NetworKit::Partition NetworKit::BinaryEdgeListPartitionReader::read(const std::string& path) {
	switch (width) {
		case 4:
			return readPartition<4>(std::vector<std::string>(1, path), firstNode);
		default:
			return readPartition<8>(std::vector<std::string>(1, path), firstNode);
	};
}

NetworKit::Partition NetworKit::BinaryEdgeListPartitionReader::read(const std::vector<std::string>& paths) {
	switch (width) {
		case 4:
			return readPartition<4>(paths, firstNode);
		default:
			return readPartition<8>(paths, firstNode);
	};
}
