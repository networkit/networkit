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

public:
	FastSNAPGraphReader() = default;

	virtual Graph read(const std::string& path) override;
};

} /* namespace NetworKit */

#endif /* FASTSNAPGRAPHREADER_H_ */
