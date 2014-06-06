/*
 * FastMETISGraphReader.h
 *
 *  Created on: 24.04.2014
 *      Author: Maximilian Vogel
 */

#ifndef FASTMETISGRAPHREADER_H_
#define FASTMETISGRAPHREADER_H_

#include "GraphReader.h"
#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"

namespace NetworKit {

class FastMETISGraphReader : public NetworKit::GraphReader {
public:
	FastMETISGraphReader() = default;

	virtual Graph read(const std::string& path) override;
};

} /* namespace NetworKit */
#endif /* FASTMETISGRAPHREADER_H_ */
