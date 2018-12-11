
/*
 * ThrillGraphBinaryReader.h
 *
 * @author Michael Hamann <michael.hamann@kit.edu>
 */

#ifndef THRILLGRAPHBINARYREADER_H_
#define THRILLGRAPHBINARYREADER_H_

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>


#include "GraphReader.h"

namespace NetworKit {

/**
 * Reads a graph format consisting of a serialized DIA of vector<uint32_t> from thrill.
 */
class ThrillGraphBinaryReader: public GraphReader {

public:

	/**
	 * When the number of nodes is given, reading the graph is more efficient.
	 *
	 * @param[in]  n  The number of nodes
	 */
	ThrillGraphBinaryReader(count n = 0);

	/**
	 * Given the path of an input file, read the graph contained.
	 *
	 * @param[in]	path	input file path
	 */
	Graph read(const std::string& path) override;
	Graph read(const std::vector<std::string>& path);
private:
	const count n;
};

} /* namespace NetworKit */
#endif /* THRILLGRAPHBINARYREADER_H_ */
