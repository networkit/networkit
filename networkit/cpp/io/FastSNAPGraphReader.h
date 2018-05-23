/*
 * FastSNAPGraphReader.h
 *
 *  Created on: 04.05.2018
 *      Author: Alexander van der Grinten
 */
#ifndef FASTSNAPGRAPHREADER_H_
#define FASTSNAPGRAPHREADER_H_

#include <unordered_map>

#include "../graph/Graph.h"
#include "GraphReader.h"

namespace NetworKit {

/**
 * @ingroup io
 */
class FastSNAPGraphReader : public NetworKit::GraphReader {
private:
	std::unordered_map<node, node> nodeIdMap;
	bool directed;
	count nodeCount;

public:
	/**
	 * Default constructor for the FastSNAPGraphReader.
	 * NOTE: Keep in mind that many SNAP graphs do not have consecutive node ids.
	 * Passing a maxNode might result in a higher memory use for the output graph
	 * but reduces runtime. Leaving this parameter on its default value, the reader
	 * is forced to remap the nodes and thereby decreases the memory use on cost of
	 * slighty higher runtime.
	 *
	 * @param[in]	directed	reads in the graph as directed, if set to true
	 * @param[in]	maxNode	maxNode is used to preallocated memory for the number of nodes
	 */
	FastSNAPGraphReader(const bool& directed = false, const count& nodeCount = 0);

	/**
	 * Given the path of an input file, read the graph contained.
	 *
	 * @param[in]	path	input file path
	 */
	virtual Graph read(const std::string& path) override;
};

} /* namespace NetworKit */

#endif /* FASTSNAPGRAPHREADER_H_ */
