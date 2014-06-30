/*
 * SNAPGraphReader.h
 *
 *  Created on: 19.05.2014
 *      Author: Maximilian Vogel
 */

#ifndef SNAPGRAPHREADER_H_
#define SNAPGRAPHREADER_H_

//#include <unordered_set>
//#include <vector>
//#include <fstream>

#include <unordered_map>

#include "../graph/Graph.h"
#include "../structures/Partition.h"
#include "GraphReader.h"

namespace NetworKit {

/**
 * @ingroup io
 */
class SNAPGraphReader : public NetworKit::GraphReader {
protected:
	std::unordered_map<node,node> mapNodeIds;

public:

	/** Default constructor */
	SNAPGraphReader() = default;

	virtual Graph read(const std::string& path) override;

	std::unordered_map<node,node> getNodeIdMap();

};

} /* namespace NetworKit */
#endif /* SNAPGRAPHREADER_H_ */
