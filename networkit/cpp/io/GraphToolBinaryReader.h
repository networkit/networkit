/*
 * GraphToolBinaryReader.h
 *
 *  Created on: 02.12.14
 *      Author: Maximilian Vogel
 */

#ifndef GRAPHTOOLBINARYREADER_H_
#define GRAPHTOOLBINARYREADER_H_

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>


#include "GraphReader.h"

namespace NetworKit {

/**
 * Reads graphs from files in the binary format defined by graph-tool[1].
 * [1]: http://graph-tool.skewed.de/static/doc/gt_format.html
 */
class GraphToolBinaryReader: public NetworKit::GraphReader {

public:

	GraphToolBinaryReader() = default; //nullary constructor for Python shell

	/**
	 * Given the path of an input file, read the graph contained.
	 *
	 * @param[in]	path	input file path
	 */
	Graph read(const std::string& path);

protected:
	bool littleEndianness; 

private:
	void addOutNeighbours(std::ifstream& file, uint64_t numNodes, Graph& G);

	std::vector<std::vector<uint64_t>> getOutNeighbours(std::ifstream& file, uint64_t numNodes);

	uint8_t getAdjacencyWidth(uint64_t n);

	uint64_t getNumNodes(std::ifstream& file);

	bool getDirected(std::ifstream& file);

	std::string readComment(std::ifstream& file);

	bool checkHeader(std::ifstream& file);

	template<typename Type>
	Type readType(std::ifstream& file, int width) {
		Type val = 0;
		uint8_t* bytes = new uint8_t[width];
		file.read((char*)bytes,width);
		if (this->littleEndianness) {
			for (int i = 0; i < width;++i) {
				val |= ((Type)bytes[i] << (i*8));
			}
		} else {
			for (int i = 0; i < width; ++i) {
				val |= ((Type)bytes[i] << (width-1-i)*8);
			}		
		}
		delete[] bytes;
		return val;
	};

};

} /* namespace NetworKit */
#endif /* GRAPHTOOLBINARYREADER_H_ */
