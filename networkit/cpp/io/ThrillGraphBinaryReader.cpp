#include "ThrillGraphBinaryReader.h"
#include "../graph/GraphBuilder.h"
#include <algorithm>

NetworKit::ThrillGraphBinaryReader::ThrillGraphBinaryReader(count n) : n(n) {};

namespace {
	template <typename stream_t>
	uint64_t GetVarint(stream_t &is) {
		auto get_byte = [&is]() -> uint8_t {
			uint8_t result;
			is.read(reinterpret_cast<char*>(&result), 1);
			return result;
		};
		uint64_t u, v = get_byte();
		if (!(v & 0x80)) return v;
		v &= 0x7F;
		u = get_byte(), v |= (u & 0x7F) << 7;
		if (!(u & 0x80)) return v;
		u = get_byte(), v |= (u & 0x7F) << 14;
		if (!(u & 0x80)) return v;
		u = get_byte(), v |= (u & 0x7F) << 21;
		if (!(u & 0x80)) return v;
		u = get_byte(), v |= (u & 0x7F) << 28;
		if (!(u & 0x80)) return v;
		u = get_byte(), v |= (u & 0x7F) << 35;
		if (!(u & 0x80)) return v;
		u = get_byte(), v |= (u & 0x7F) << 42;
		if (!(u & 0x80)) return v;
		u = get_byte(), v |= (u & 0x7F) << 49;
		if (!(u & 0x80)) return v;
		u = get_byte(), v |= (u & 0x7F) << 56;
		if (!(u & 0x80)) return v;
		u = get_byte();
		if (u & 0xFE)
			throw std::overflow_error("Overflow during varint64 decoding.");
		v |= (u & 0x7F) << 63;
		return v;
	}
};


NetworKit::Graph NetworKit::ThrillGraphBinaryReader::read(const std::string &path) {
	return read(std::vector<std::string>(1, path));
};

NetworKit::Graph NetworKit::ThrillGraphBinaryReader::read(const std::vector<std::string> &paths) {
	GraphBuilder gb(n);

	if (!paths.empty()) {
		std::ifstream is;
		count file_index = 0;

		auto next_input = [&]() {
			is.close();
			if (file_index < paths.size()) {
				is.open(paths[file_index++]);
			}
		};

		next_input();

		node max_id = 0;

		for (node u = 0; is.good() && is.is_open(); ++u) {
			// Add node if it does not exist yet, only one may be missing
			if (u >= gb.upperNodeIdBound()) {
				gb.addNode();
			}

			for (count deg = GetVarint(is); deg > 0 && is.good(); --deg) {
				uint32_t v;
				if (!is.read(reinterpret_cast<char*>(&v), 4)) {
					throw std::runtime_error("I/O error while reading next neighbor");
				}

				max_id = std::max<node>(max_id, v);

				gb.addHalfEdge(u, v);
			}

			if (is.is_open() && (is.peek() == std::char_traits<char>::eof() || !is.good())) {
				next_input();
			}
		}

		if (max_id >= gb.upperNodeIdBound()) {
			throw std::runtime_error("Maximum read node id larger than number of nodes read.");
		}
	}

	return gb.toGraph(true, true);
};
