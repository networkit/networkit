/*
 * GraphToolBinaryWriter.h
 *
 *  Created on: 02.12.14
 *      Author: Maximilian Vogel
 */

#ifndef GRAPHTOOLBINARYWRITER_H_
#define GRAPHTOOLBINARYWRITER_H_

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>


#include "GraphWriter.h"

namespace NetworKit {

/**
 *
 */
class GraphToolBinaryWriter: public NetworKit::GraphWriter {

public:

	//GraphToolBinaryWriter() = default; //nullary constructor for Python shell

	GraphToolBinaryWriter(bool littleEndianness = true);

	/**
	 * Given the path of an input file, read the graph contained.
	 *
	 * @param[in]	path	input file path
	 */
	void write(Graph& G, const std::string& path);

protected:
	bool littleEndianness;


private:
	void writeAdjacencies(std::ofstream& file, Graph& G);

	uint8_t getAdjacencyWidth(uint64_t n);

	//uint64_t getNumNodes(std::istream& file);

	//bool getDirected(std::istream& file);

	void writeComment(std::ofstream& file);

	void writeHeader(std::ofstream& file);

	template<typename Type>
	void writeType(std::ofstream& file, int width, Type val) {
		//DEBUG("writing ",val, "with width ", width, " to file");
		uint8_t* bytes = new uint8_t[width];
		if (this->littleEndianness) {
			for (int i = 0; i < width;++i) {
				bytes[i] = (val >> (i*8)) & 0xFF;
			}
		} else {
			for (int i = 0; i < width; ++i) {
				bytes[i] = (val >> ((width-1-i)*8)) & 0xFF;
			}		
		}
		file.write((char*)bytes,width);
		delete[] bytes;
	};

};

} /* namespace NetworKit */
#endif /* GRAPHTOOLBINARYWRITER_H_ */
