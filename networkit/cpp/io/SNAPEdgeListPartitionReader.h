/*
 * SNAPEdgeListPartitionReader.h
 *
 *  Created on: Jun 20, 2013
 *      Author: forigem
 */

#ifndef SNAPEDGELISTPARTITIONREADER_H_
#define SNAPEDGELISTPARTITIONREADER_H_

#include <set>
#include <vector>
#include <fstream>
#include <unordered_map>

#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"
#include "../structures/Cover.h"

namespace NetworKit {

/**
 * @ingroup io
 * Reads the clustering files from the SNAP collection.
 */
class SNAPEdgeListPartitionReader {
public:
	virtual Cover read(std::string path, std::unordered_map<node,node>& mapNodeIds, Graph& G);

//	virtual Partition readWithInfo(std::string path, count nNodes);

};

} /* namespace NetworKit */
#endif /* SNAPEDGELISTPARTITIONREADER_H_ */
