/*
 * NetworkitBinaryWriter.cpp
 *
 *@author Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#include <networkit/io/NetworkitBinaryWriter.hpp>
#include <networkit/auxiliary/Enforce.hpp>
#include <tlx/math/clz.hpp>
#include <fstream>
#include <cstring>

namespace NetworKit {

NetworkitBinaryWriter::NetworkitBinaryWriter(uint64_t chunks, NetworkitBinaryWeights weightsType): chunks(chunks), weightsType(weightsType) {}

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
	nkbg::WEIGHT_FORMAT weightFormat;

	auto detectWeightsType = [&] () {
		weightFormat = nkbg::WEIGHT_FORMAT::unsignedFormat;
		bool isUnsigned = true;
		bool fitsIntoInt64 = true;
		bool fitsIntoFloat = true;
		G.forNodes([&](node n) {
			G.forNeighborsOf(n,[&](node v, edgeweight w) {
				int64_t weight = w*10;
				if(w < 0)
					isUnsigned = false;
				if(w != static_cast<int64_t>(w))
					fitsIntoInt64 = false;
				if(w != static_cast<float>(w))
					fitsIntoFloat = false;
			});
		});
		if(fitsIntoInt64) {
			if(isUnsigned) {
				weightFormat = nkbg::WEIGHT_FORMAT::unsignedFormat;
			} else {
				weightFormat = nkbg::WEIGHT_FORMAT::signedFormat;
			}
		} else {
			if(fitsIntoFloat) {
				weightFormat = nkbg::WEIGHT_FORMAT::floatFormat;
			} else {
				weightFormat = nkbg::WEIGHT_FORMAT::doubleFormat;
			}
		}
	};

	switch(weightsType) {
		case NetworkitBinaryWeights::none:
			weightFormat = nkbg::WEIGHT_FORMAT::none;
		break;
		case NetworkitBinaryWeights::autoDetect:
			detectWeightsType();
		break;
		case NetworkitBinaryWeights::unsignedFormat:
			weightFormat = nkbg::WEIGHT_FORMAT::unsignedFormat;
		break;
		case NetworkitBinaryWeights::signedFormat:
			weightFormat = nkbg::WEIGHT_FORMAT::signedFormat;
		break;
		case NetworkitBinaryWeights::doubleFormat:
			weightFormat = nkbg::WEIGHT_FORMAT::doubleFormat;
		break;
		case NetworkitBinaryWeights::floatFormat:
			weightFormat = nkbg::WEIGHT_FORMAT::floatFormat;
		break;
	}

	nkbg::Header header = {};
	uint64_t adjListSize = 0;
	uint64_t adjTransposeSize = 0;

	auto setFeatures = [&] () {
		header.features |= (G.isDirected() & nkbg::DIR_MASK);
		header.features |= ((static_cast<uint64_t>(weightFormat) << nkbg::WGHT_SHIFT) & nkbg::WGHT_MASK);
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
		outfile.write(reinterpret_cast<char*>(&header.offsetWeightLists), sizeof(uint64_t));
		outfile.write(reinterpret_cast<char*>(&header.offsetWeightTranspose), sizeof(uint64_t));
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
	uint64_t adjWeightSize = 0;
	uint64_t transpWeightSize = 0;
	std::vector<uint64_t> nrOutNbrs;
	std::vector<uint64_t> nrInNbrs;
	std::vector<size_t> adjOffsets;	//Prefix sum of size encoded adj arrays
	std::vector<size_t> transpOffsets;	//Prefix sum of encoded transposed adj arrays
	std::vector<size_t> adjWghtOffsets;	//Prefix sum of size encoded adj weights
	std::vector<size_t> transpWghtOffsets;	//Prefix sum of encoded transposed weights

	for(uint64_t c = 0; c < chunks; c++) {
		for(uint64_t n = firstInChunk[c]; n < firstInChunk[c+1]; n++) {
			uint64_t outNbrs = 0;
			uint64_t inNbrs = 0;
			uint8_t tmp [10];
			uint64_t weight;
			if(!G.isDirected()){
				G.forNeighborsOf(n,[&](node v, edgeweight w) {
					if(v <= n) {
						outNbrs++;
						adjSize += encode(v, tmp);
						if(weightFormat == nkbg::WEIGHT_FORMAT::unsignedFormat) {
							adjWeightSize += encode(w,tmp);
						} else if (weightFormat == nkbg::WEIGHT_FORMAT::doubleFormat) {
							adjWeightSize += sizeof(double);
						} else if(weightFormat == nkbg::WEIGHT_FORMAT::signedFormat) {
							uint64_t weight = encodeZigzag(w);
							adjWeightSize += encode(weight,tmp);
						} else if(weightFormat == nkbg::WEIGHT_FORMAT::floatFormat) {
							adjWeightSize += sizeof(float);
						}
					} else if (v >= n) {
						inNbrs++;
						transpSize += encode(v, tmp);
						if(weightFormat == nkbg::WEIGHT_FORMAT::unsignedFormat) {
							transpWeightSize += encode(w,tmp);
						} else if(weightFormat == nkbg::WEIGHT_FORMAT::doubleFormat) {
							transpWeightSize += sizeof(double);
						} else if(weightFormat == nkbg::WEIGHT_FORMAT::signedFormat) {
							uint64_t weight = encodeZigzag(w);
							transpWeightSize += encode(weight,tmp);
						} else if(weightFormat == nkbg::WEIGHT_FORMAT::floatFormat) {
							transpWeightSize += sizeof(float);
						}
					}
				});
			} else {
				G.forNeighborsOf(n, [&] (node v, edgeweight w) {
					outNbrs++;
					adjSize += encode(v, tmp);
					if(weightFormat == nkbg::WEIGHT_FORMAT::unsignedFormat) {
						adjWeightSize += encode(w,tmp);
					} else if(weightFormat == nkbg::WEIGHT_FORMAT::doubleFormat) {
						adjWeightSize += sizeof(double);
					} else if(weightFormat == nkbg::WEIGHT_FORMAT::signedFormat) {
						uint64_t weight = encodeZigzag(w);
						adjWeightSize += encode(weight,tmp);
					} else if(weightFormat == nkbg::WEIGHT_FORMAT::floatFormat) {
						adjWeightSize += sizeof(float);
					}
				});
				G.forInNeighborsOf(n, [&] (node v, edgeweight w) {
					inNbrs++;
					transpSize += encode(v, tmp);
					if(weightFormat == nkbg::WEIGHT_FORMAT::unsignedFormat) {
						transpWeightSize += encode(w,tmp);
					} else if(weightFormat == nkbg::WEIGHT_FORMAT::doubleFormat) {
						transpWeightSize += sizeof(double);
					} else if(weightFormat == nkbg::WEIGHT_FORMAT::signedFormat) {
						uint64_t weight = encodeZigzag(w);
						transpWeightSize += encode(weight,tmp);
					} else if(weightFormat == nkbg::WEIGHT_FORMAT::floatFormat) {
						transpWeightSize += sizeof(float);
					}
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
		adjWghtOffsets.push_back(adjWeightSize);
		transpWghtOffsets.push_back(transpWeightSize);
	}
	// Write header.
	strncpy(header.magic,"nkbg002",8);
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
	if(weightFormat != nkbg::WEIGHT_FORMAT::none) {
			header.offsetWeightLists = header.offsetAdjTranspose
			+ (chunks - 1) * sizeof(uint64_t) //transpOffsets
			+ sizeof(uint64_t) //adjTranspSize
			+ transpOffsets.back(); //Size of data
		header.offsetWeightTranspose = header.offsetWeightLists
			+ (chunks - 1) * sizeof(uint64_t)
			+ adjWghtOffsets.back();
	} else {
		header.offsetWeightLists = 0;
		header.offsetWeightTranspose = 0;
	}
	writeHeader();
	// Write base data.
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
			uint64_t weightSize;
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
	// Write adj weights.
	for (uint64_t c = 1; c < chunks; c++) {
		outfile.write(reinterpret_cast<char*>(&adjWghtOffsets[c-1]), sizeof(uint64_t));
	}
	G.forNodes([&](node u) {
		uint8_t tmp [10];
		G.forNeighborsOf(u,[&](node v,edgeweight w){
			uint64_t weightSize;
			if(!G.isDirected()) {
				if(v <= u) {
					switch(weightFormat) {
						case nkbg::WEIGHT_FORMAT::unsignedFormat:
							weightSize = encode(w, tmp);
							outfile.write(reinterpret_cast<char*>(tmp), weightSize);
						break;
						case nkbg::WEIGHT_FORMAT::doubleFormat:
							{
								outfile.write(reinterpret_cast<char*>(&w), sizeof(double));
							}
						break;
						case nkbg::WEIGHT_FORMAT::signedFormat:
							{
								uint64_t weight = encodeZigzag(w);
								weightSize = encode(weight,tmp);
								outfile.write(reinterpret_cast<char*>(tmp), weightSize);
							}
						break;
						case nkbg::WEIGHT_FORMAT::floatFormat:
							{
								float weight = w;
								outfile.write(reinterpret_cast<char*>(&weight), sizeof(float));
							}
						break;
					}
				}
			} else {
				switch(weightFormat) {
					case nkbg::WEIGHT_FORMAT::unsignedFormat:
						{
							weightSize = encode(w, tmp);
							outfile.write(reinterpret_cast<char*>(tmp), weightSize);
						}
					break;
					case nkbg::WEIGHT_FORMAT::doubleFormat:
						{
							outfile.write(reinterpret_cast<char*>(&w), sizeof(double));
						}
					break;
					case nkbg::WEIGHT_FORMAT::signedFormat:
						{
							int64_t weight = encodeZigzag(w);
							weightSize = encode(weight,tmp);
							outfile.write(reinterpret_cast<char*>(tmp), weightSize);
						}
					break;
					case nkbg::WEIGHT_FORMAT::floatFormat:
						{
							float weight = w;
							outfile.write(reinterpret_cast<char*>(&weight), sizeof(float));
						}
					break;
				}
			}
		});
	});

	// Write transpose weights.
	for (uint64_t c = 1; c < chunks; c++) {
		outfile.write(reinterpret_cast<char*>(&transpWghtOffsets[c-1]), sizeof(uint64_t));
	}
	G.forNodes([&](node u) {
		uint8_t tmp [10];
		uint64_t weightSize;
		if(!G.isDirected()) {
			G.forNeighborsOf(u,[&](node v, edgeweight w){
				if(v >= u) {
					switch(weightFormat) {
						case nkbg::WEIGHT_FORMAT::unsignedFormat:
							{
								weightSize = encode(w, tmp);
								outfile.write(reinterpret_cast<char*>(tmp), weightSize);
							}
						break;
						case nkbg::WEIGHT_FORMAT::doubleFormat:
							{
								outfile.write(reinterpret_cast<char*>(&w), sizeof(double));
							}
						break;
						case nkbg::WEIGHT_FORMAT::signedFormat:
							{
								uint64_t weight = encodeZigzag(w);
								weightSize = encode(weight, tmp);
								outfile.write(reinterpret_cast<char*>(tmp), weightSize);
							}
						break;
						case nkbg::WEIGHT_FORMAT::floatFormat:
							{
								float weight = w;
								outfile.write(reinterpret_cast<char*>(&weight), sizeof(float));
							}
						break;
					}
				}
			});
		} else {
			G.forInNeighborsOf(u,[&](node v, edgeweight w){
				switch(weightFormat) {
					case nkbg::WEIGHT_FORMAT::unsignedFormat:
						{
							weightSize = encode(w, tmp);
							outfile.write(reinterpret_cast<char*>(tmp), weightSize);
						}
					break;
					case nkbg::WEIGHT_FORMAT::doubleFormat:
						{
							outfile.write(reinterpret_cast<char*>(&w), sizeof(double));
						}
					break;
					case nkbg::WEIGHT_FORMAT::signedFormat:
						{
							uint64_t weight = encodeZigzag(w);
							weightSize = encode(weight, tmp);
							outfile.write(reinterpret_cast<char*>(tmp), weightSize);
						}
					break;
					case nkbg::WEIGHT_FORMAT::floatFormat:
						{
							float weight = w;
							outfile.write(reinterpret_cast<char*>(&weight), sizeof(float));
						}
					break;
				}
			});
		}
	});
	INFO("Written graph to ", path);
}
} /*namespace */
