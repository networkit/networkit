#include "BinaryEdgeListPartitionWriter.h"
#include <fstream>

NetworKit::BinaryEdgeListPartitionWriter::BinaryEdgeListPartitionWriter(NetworKit::node firstNode, uint8_t width) : firstNode(firstNode), width(width) {
	if (width != 4 && width != 8) {
 		throw std::runtime_error("Width must be 4 or 8");
	}
}

namespace {
	template <uint8_t width>
	void writePartition(NetworKit::Partition &zeta, const std::string &path) {
		using write_t = typename std::conditional<width == 4, uint32_t, uint64_t>::type;
		static_assert(sizeof(write_t) == width, "write_t does not have the expected width");
		
		std::ofstream os(path, std::ios::trunc | std::ios::binary);
		
		zeta.forEntries([&](NetworKit::index u, NetworKit::index p) {
				write_t uw = u;
				write_t pw = p;
				os.write(reinterpret_cast<const char*>(&uw), width);
				os.write(reinterpret_cast<const char*>(&pw), width);
		    });
	}
}

void NetworKit::BinaryEdgeListPartitionWriter::write( NetworKit::Partition &zeta, const std::string &path ) const {
	switch (width) {
	case 4:
		writePartition<4>(zeta, path);
		break;
	default:
		writePartition<8>(zeta, path);
		break;
	}
}
