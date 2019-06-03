/*
 * NetworkitBinaryReader.cpp
 *
 *@author Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#include <networkit/io/NetworkitBinaryReader.hpp>

namespace NetworKit {

size_t NetworkitBinaryReader::decode(const uint8_t* data, uint64_t& value) {
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

int64_t NetworkitBinaryReader::decodeZigzag(uint64_t value) {
	 return (value >> 1) ^ (-(value & 1));

}

Graph NetworkitBinaryReader::read(const std::string& path) {
	Header header = {};
	uint64_t weightFormat = 0;

	MemoryMappedFile mmfile(path);
	auto it = mmfile.cbegin();
	auto end = mmfile.cend();

	auto readHeader = [&] () {
		memcpy(&header.magic, it, sizeof(uint64_t));
		it += sizeof(uint64_t);
		memcpy(&header.checksum, it, sizeof(uint64_t));
		it += sizeof(uint64_t);
		memcpy(&header.features, it, sizeof(uint64_t));
		it += sizeof(uint64_t);
		memcpy(&header.nodes, it, sizeof(uint64_t));
		it += sizeof(uint64_t);
		memcpy(&header.chunks, it, sizeof(uint64_t));
		it += sizeof(uint64_t);
		memcpy(&header.offsetBaseData, it, sizeof(uint64_t));
		it += sizeof(uint64_t);
		memcpy(&header.offsetAdjLists, it, sizeof(uint64_t));
		it += sizeof(uint64_t);
		memcpy(&header.offsetWeights, it, sizeof(uint64_t));
		it += sizeof(uint64_t);
	};

	auto checkHeader = [&] () {
		if(memcmp("nkbg000", header.magic, 8)) {
			throw std::runtime_error("Reader expected another magic value");
		} else {
			directed = (header.features & DIR_MASK);
			weightFormat = (header.features & WGHT_MASK) >> WGHT_SHIFT;
		}
	};

	readHeader();
	checkHeader();

	nodes = header.nodes;
	DEBUG("# nodes here = ", nodes);
	chunks = header.chunks;
	DEBUG("# chunks here = ", chunks);
	weighted = header.offsetWeights;
	Graph G(nodes, weighted, directed);

	// Read base data.
	std::vector<uint64_t> prefSum;
	const char *baseIt = mmfile.cbegin() + header.offsetBaseData;
	for(uint64_t i = 0; i < nodes; i++) {
		uint64_t sum;
		memcpy(&sum, baseIt, sizeof(uint64_t));
		baseIt += sizeof(uint64_t);
		uint64_t deleted = (sum & DELETED_BIT);
		if(!deleted) {
			sum &= SIZE_MASK;
			prefSum.push_back(sum);
		} else {
			G.removeNode(i);
		}
	}

	std::vector<uint64_t> firstVert;
	firstVert.push_back(0);
	for(uint64_t ch = 1; ch < chunks; ch++) {
		uint64_t first;
		memcpy(&first, baseIt, sizeof(uint64_t));
		baseIt += sizeof(uint64_t);
		firstVert.push_back(first);
	}
	firstVert.push_back(nodes);

	// Read adjacency lists.
	const char *adjIt = mmfile.cbegin() + header.offsetAdjLists;

	auto constructGraph = [&] (uint64_t c) {
		node vertex = firstVert[c];
		uint64_t off = 0;
		if(vertex) {
			memcpy(&off, adjIt + (c-1)* sizeof(uint64_t), sizeof(uint64_t));
		}
		off += (chunks - 1) * sizeof(uint64_t);
		uint64_t n = firstVert[c+1] - firstVert[c];

		for (uint64_t i = 0; i < n; i++) {
			uint64_t nbrs;
			uint64_t curr = vertex+i;
			if(curr > 0) {
				nbrs = prefSum[curr] - prefSum[curr-1];
			} else {
				nbrs = prefSum[0];
			}
			if(!directed) {
				G.preallocateUndirected(curr, nbrs);
			}
			for (uint64_t j = 0; j < nbrs; j++) {
				uint64_t add;
				size_t rd = NetworkitBinaryReader::decode(reinterpret_cast<const uint8_t*>(adjIt + off), add);
				if(!G.isDirected()) {
					if(add <= curr) {
						G.addEdge(curr, add);
					}
				} else {
					G.addEdge(curr,add);
				}
				off+=rd;
			}
		}
	};

	// create graph
	for(uint64_t c = 0; c < chunks; c++) {
		constructGraph(c);
	}
	return G;
}
} /* namespace */


