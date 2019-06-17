/*
 * NetworkitBinaryWriter.cpp
 *
 *@author Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#include <networkit/io/NetworkitBinaryWriter.hpp>

namespace NetworKit {

NetworkitBinaryWriter::NetworkitBinaryWriter(uint64_t chunks): chunks(chunks) {}

size_t NetworkitBinaryWriter::encode(uint64_t value, uint8_t* buffer) {
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

uint64_t NetworkitBinaryWriter::encodeZigzag(int64_t value) {
	return (value << 1) ^ (value >> 31);
}

void NetworkitBinaryWriter::write(const Graph &G, const std::string& path) {

	std::ofstream outfile(path, std::ios::binary);
	Aux::enforceOpened(outfile);

	uint64_t weightFormat = 0;
	Header header = {};

	auto setFeatures = [&] () {
		header.features |= (G.isDirected() & DIR_MASK);
		header.features |= ((weightFormat << WGHT_SHIFT) & WGHT_MASK);
	};

	auto writeHeader = [&] () {
		outfile.write(header.magic, 8);
		outfile.write(reinterpret_cast<char*>(&header.checksum), sizeof(uint64_t));
		outfile.write(reinterpret_cast<char*>(&header.features), sizeof(uint64_t));
		outfile.write(reinterpret_cast<char*>(&header.nodes), sizeof(uint64_t));
		outfile.write(reinterpret_cast<char*>(&header.chunks), sizeof(uint64_t));
		outfile.write(reinterpret_cast<char*>(&header.offsetBaseData), sizeof(uint64_t));
		outfile.write(reinterpret_cast<char*>(&header.offsetAdjLists), sizeof(uint64_t));
		outfile.write(reinterpret_cast<char*>(&header.offsetAdjTranspose), sizeof(uint64_t));
		outfile.write(reinterpret_cast<char*>(&header.offsetWeights), sizeof(uint64_t));
	};

	nodes = G.numberOfNodes();
	if (nodes < chunks) {
		chunks = nodes;
		INFO("reducing chunks to ", chunks, " chunks");
	}

	// Compute first vertex of each chunk.
	std::vector<uint64_t> firstInChunk;
	firstInChunk.push_back(0);
	uint64_t firstNode = 0;
	for(uint64_t c = 1; c < chunks; c++) {
		firstNode += (nodes/chunks);
		firstInChunk.push_back(firstNode);
	}
	firstInChunk.push_back(nodes);

	// Compute encoded size of arrays and store in vector.
	uint64_t adjSize = 0; 
	std::vector<size_t> adjOffsets;	//  pref sum of size encoded adj arrays
	for(uint64_t c = 0; c < chunks; c++) {
		for(uint64_t n = firstInChunk[c]; n < firstInChunk[c+1]; n++) {
			G.forNeighborsOf(n,[&](node v) {
				uint8_t tmp [10];
				adjSize += encode(v, tmp);
			});
		}
		adjOffsets.push_back(adjSize);
	}

	// Write header.
	strncpy(header.magic,"nkbg000",8);
	header.checksum = 0;
	setFeatures();
	header.nodes = nodes;
	header.chunks = chunks;
	header.offsetBaseData = sizeof(Header);
	header.offsetAdjLists = header.offsetBaseData
			+ nodes * sizeof(uint64_t) // prefixSum.
			+ (chunks - 1) * sizeof(uint64_t); // firstVertex.
	header.offsetAdjTranspose = 0;
	header.offsetWeights = 0;

	writeHeader();

	// Write base data.
	uint64_t sum = 0;
	G.forNodes([&] (node u) {
		sum += G.degree(u);
		assert(!(sum & ~SIZE_MASK));
		if(!G.hasNode(u)) {
			sum |= DELETED_BIT;
		}
		outfile.write(reinterpret_cast<char*>(&sum), sizeof(uint64_t));
	});

	assert(!firstInChunk[0]);
	for (uint64_t c = 1; c < chunks; c++) {
		outfile.write(reinterpret_cast<char*>(&firstInChunk[c]), sizeof(uint64_t));
	}

	// Write adjacency data.
	for (uint64_t c = 1; c < chunks; c++) {
		outfile.write(reinterpret_cast<char*>(&adjOffsets[c-1]), sizeof(uint64_t));
	}

	G.forNodes([&](node u) {
		G.forNeighborsOf(u,[&](node v){
			uint8_t tmp [10];
			uint64_t nodeSize = encode(v, tmp);
			outfile.write(reinterpret_cast<char*>(tmp), nodeSize);
		});
	});

	INFO("Written graph to ", path);
}
} /*namespace */
