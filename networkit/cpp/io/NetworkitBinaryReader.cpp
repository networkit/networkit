/*
 * NetworkitBinaryReader.cpp
 *
 *@author Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#include <networkit/io/NetworkitBinaryReader.hpp>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>
#include <networkit/io/MemoryMappedFile.hpp>
#include <tlx/math/ffs.hpp>
#include <vector>
#include <fstream>
#include <string.h>
#include <atomic>

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
	nkbg::Header header = {};
	uint64_t weightFormat = 0;

	MemoryMappedFile mmfile(path);
	auto it = mmfile.cbegin();

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
		memcpy(&header.offsetAdjTranspose, it, sizeof(uint64_t));
		it += sizeof(uint64_t);
		memcpy(&header.offsetWeights, it, sizeof(uint64_t));
		it += sizeof(uint64_t);
	};

	auto checkHeader = [&] () {
		if(memcmp("nkbg001", header.magic, 8)) {
			throw std::runtime_error("Reader expected another magic value");
		} else {
			directed = (header.features & nkbg::DIR_MASK);
			weightFormat = (header.features & nkbg::WGHT_MASK) >> nkbg::WGHT_SHIFT;
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
	std::vector<uint8_t> nodeFlags;
	const char *baseIt = mmfile.cbegin() + header.offsetBaseData;
	for(uint64_t i = 0; i < nodes; i++) {
		uint8_t flag;
		memcpy(&flag, baseIt, sizeof(uint8_t));
		baseIt += sizeof(uint8_t);
		if (!(flag & nkbg::DELETED_BIT)) {
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
	const char *transpIt = mmfile.cbegin() + header.offsetAdjTranspose;
	uint64_t adjListSize;
	memcpy(&adjListSize, adjIt + (chunks -1) * sizeof(uint64_t), sizeof(uint64_t));
	uint64_t transposeListSize;
	memcpy(&transposeListSize, transpIt + (chunks -1) * sizeof(uint64_t), sizeof(uint64_t));

	if(!directed) {
		assert(adjListSize == transposeListSize);
	} 
	G.setEdgeCount(unsafe, adjListSize);

	std::atomic<count> selfLoops{0};
	auto constructGraph = [&] (uint64_t c) {
		node vertex = firstVert[c];
		uint64_t off = 0;
		uint64_t transpOff = 0;
		if(vertex) {
			memcpy(&off, adjIt + (c-1)* sizeof(uint64_t), sizeof(uint64_t));
			memcpy(&transpOff, transpIt + (c-1)* sizeof(uint64_t), sizeof(uint64_t));
		}
		off += (chunks - 1) * sizeof(uint64_t);
		transpOff += (chunks - 1) * sizeof(uint64_t);
		off += sizeof(uint64_t);
		transpOff += sizeof(uint64_t);
		uint64_t n = firstVert[c+1] - firstVert[c];

		for (uint64_t i = 0; i < n; i++) {
			uint64_t curr = vertex+i;
			uint64_t outNbrs;
			off += NetworkitBinaryReader::decode(reinterpret_cast<const uint8_t*>(adjIt + off), outNbrs);
			uint64_t inNbrs;
			transpOff += NetworkitBinaryReader::decode(reinterpret_cast<const uint8_t*>(transpIt + transpOff), inNbrs);
			if(!directed) {
				G.preallocateUndirected(curr, outNbrs+inNbrs);
			} else  {
				G.preallocateDirected(curr, outNbrs, inNbrs);
			}
			//Read adjacency lists.
			for (uint64_t j = 0; j < outNbrs; j++) {
				uint64_t add;
				off += NetworkitBinaryReader::decode(reinterpret_cast<const uint8_t*>(adjIt + off), add);
				if(!directed) {
					G.addPartialEdge(unsafe, curr, add);
				} else {
					G.addPartialOutEdge(unsafe, curr, add);
				}
				if(curr == add) {
					selfLoops++;
				}
			}
			//Read transpose lists.
			for (uint64_t j = 0; j < inNbrs; j++) {
				uint64_t add;
				transpOff += NetworkitBinaryReader::decode(reinterpret_cast<const uint8_t*>(transpIt + transpOff), add);
				if(!directed) {
					G.addPartialEdge(unsafe, curr, add);
				} else {
					G.addPartialInEdge(unsafe, curr, add);
				}
				if(curr == add) {
					selfLoops.fetch_add(1, std::memory_order_relaxed);
				}
			}
		}
	};

	// create graph
	#pragma omp parallel for
	for(omp_index c = 0; c < static_cast<omp_index>(chunks); c++) {
		constructGraph(c);
	}
	G.setNumberOfSelfLoops(unsafe, selfLoops);
	return G;
}
} /* namespace */


