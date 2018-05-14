/*
 * FastKONECTGraphReader.h
 *
 *  Created on: 11.05.2018
 *      Author: Roman Bange
 */

#ifndef FASTKONECTGRAPHREADER_H_
#define FASTKONECTGRAPHREADER_H_

#include <unordered_map>

#include "../graph/Graph.h"
#include "GraphReader.h"

namespace NetworKit {
  class FastKONECTGraphReader : public NetworKit::GraphReader{

    public:
    	FastKONECTGraphReader() = default;

    	virtual Graph read(const std::string& path) override;
  };
} /* namespace NetworKit */
#endif /* FASTKONECTGRAPHREADER_H_ */
