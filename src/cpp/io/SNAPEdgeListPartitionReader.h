/*
 * SNAPEdgeListPartitionReader.h
 *
 *  Created on: Jun 20, 2013
 *      Author: forigem
 */

#ifndef SNAPEDGELISTPARTITIONREADER_H_
#define SNAPEDGELISTPARTITIONREADER_H_

#include <unordered_set>
#include <vector>
#include <fstream>

#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"

namespace NetworKit {

/**
 * Reads the clustering files from the SNAP collection.
 */
class SNAPEdgeListPartitionReader {
public:
	SNAPEdgeListPartitionReader();
	virtual ~SNAPEdgeListPartitionReader();

	virtual std::vector<std::set<node>>  read(std::string path);

};

} /* namespace NetworKit */
#endif /* SNAPEDGELISTPARTITIONREADER_H_ */
