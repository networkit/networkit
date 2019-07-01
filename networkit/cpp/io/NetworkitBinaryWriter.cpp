/*
 * NetworkitBinaryWriter.cpp
 *
 *@author Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#include <networkit/io/NetworkitBinaryWriter.hpp>
#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>
#include <tlx/math/clz.hpp>
#include <fstream>
#include <cstring>

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
	nkbg::Header header = {};
	uint64_t adjListSize = 0;
	uint64_t adjTransposeSize = 0;

	auto setFeatures = [&] () {
		header.features |= (G.isDirected() & nkbg::DIR_MASK);
		header.features |= ((weightFormat << nkbg::WGHT_SHIFT) & nkbg::WGHT_MASK);
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
	uint64_t transpSize = 0;
	std::vector<uint64_t> nrOutNbrs;
	std::vector<uint64_t> nrInNbrs;
	std::vector<size_t> adjOffsets;	//Prefix sum of size encoded adj arrays
	std::vector<size_t> transpOffsets;	//Prefix sum of encoded transposed adj arrays
	for(uint64_t c = 0; c < chunks; c++) {
		for(uint64_t n = firstInChunk[c]; n < firstInChunk[c+1]; n++) {
			uint64_t outNbrs = 0;
			uint64_t inNbrs = 0;
			uint8_t tmp [10];
			if(!G.isDirected()){
				G.forNeighborsOf(n,[&](node v) {
					if(v <= n) { 
						outNbrs++;
						adjSize += encode(v, tmp);
					} else if (v >= n) {
						inNbrs++;
						transpSize += encode(v, tmp);
					}
				});
			} else {
				G.forNeighborsOf(n, [&] (node v) {
					outNbrs++;
					adjSize += encode(v, tmp);
				});
				G.forInNeighborsOf(n, [&] (node v) {
					inNbrs++;
					transpSize += encode(v, tmp);
				});
			}
			adjListSize += outNbrs;
			adjSize += encode(outNbrs, tmp);
			nrOutNbrs.push_back(outNbrs);

			adjTransposeSize += inNbrs;
			transpSize += encode(inNbrs, tmp);
			nrInNbrs.push_back(inNbrs);
		}
		adjOffsets.push_back(adjSize);
		transpOffsets.push_back(transpSize);
	}

	// Write header.
	strncpy(header.magic,"nkbg001",8);
	header.checksum = 0;
	setFeatures();
	header.nodes = nodes;
	header.chunks = chunks;
	header.offsetBaseData = sizeof(nkbg::Header);
	header.offsetAdjLists = header.offsetBaseData
			+ nodes * sizeof(uint8_t) // nodeFlags.
			+ (chunks - 1) * sizeof(uint64_t); // firstVertex.
	header.offsetAdjTranspose = header.offsetAdjLists
			+ (chunks -1) * sizeof(uint64_t) // adjOffsets
			+ sizeof(uint64_t) // adjListSize
			+ adjOffsets.back(); // Size of data
	header.offsetWeights = 0;

	writeHeader();

	// Write base data.
	uint64_t sum = 0;
	G.forNodes([&] (node u) {
		uint8_t nodeFlag = 0;
		if (G.hasNode(u)) {
			nodeFlag |= nkbg::DELETED_BIT;
		}
		outfile.write(reinterpret_cast<char*>(&nodeFlag), sizeof(uint8_t));
	});

	assert(!firstInChunk[0]);
	for (uint64_t c = 1; c < chunks; c++) {
		outfile.write(reinterpret_cast<char*>(&firstInChunk[c]), sizeof(uint64_t));
	}

	// Write adjacency data.
	for (uint64_t c = 1; c < chunks; c++) {
		outfile.write(reinterpret_cast<char*>(&adjOffsets[c-1]), sizeof(uint64_t));
	}
	// Write size of list
	outfile.write(reinterpret_cast<char*>(&adjListSize), sizeof(uint64_t));
	G.forNodes([&](node u) {
		uint8_t tmp [10];
		uint64_t nbrsSize = encode(nrOutNbrs[u], tmp);
		outfile.write(reinterpret_cast<char*>(tmp), nbrsSize);
		G.forNeighborsOf(u,[&](node v){
			uint64_t nodeSize;
			if(!G.isDirected()) {
				if(v <= u) {
					nodeSize = encode(v, tmp);
					outfile.write(reinterpret_cast<char*>(tmp), nodeSize);
				}
			} else {
				nodeSize = encode(v, tmp);
				outfile.write(reinterpret_cast<char*>(tmp), nodeSize);
			}
		});
	});

	// Write transpose data.
	for (uint64_t c = 1; c < chunks; c++) {
		outfile.write(reinterpret_cast<char*>(&transpOffsets[c-1]), sizeof(uint64_t));
	}
	// Write size of transpose list.
	outfile.write(reinterpret_cast<char*>(&adjTransposeSize), sizeof(uint64_t));
	G.forNodes([&](node u) {
		uint8_t tmp [10];
		uint64_t nbrsSize = encode(nrInNbrs[u], tmp);
		outfile.write(reinterpret_cast<char*>(tmp), nbrsSize);
		uint64_t nodeSize;
		if(!G.isDirected()) {
			G.forNeighborsOf(u,[&](node v){
				if(v >= u) {
					nodeSize = encode(v, tmp);
					outfile.write(reinterpret_cast<char*>(tmp), nodeSize);
				}
			});
		} else {
			G.forInNeighborsOf(u,[&](node v){
				nodeSize = encode(v, tmp);
				outfile.write(reinterpret_cast<char*>(tmp), nodeSize);
			});
		}
	});
	INFO("Written graph to ", path);
}
} /*namespace */
