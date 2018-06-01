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
		/*
		* If the input graph has multiple edges, you can specify on how this edges are handled.
		* See GraphReader.h for a closer description of the paramters.
		*/
		FastKONECTGraphReader(MultipleEdgesHandling handlingmethod = DISCARD_EDGES);

		/**
		 * Given the path of an input file, read the graph contained.
		 *
		 * @param[in]	path	input file path
		 */
		virtual Graph read(const std::string& path) override;

	protected:
		MultipleEdgesHandling multipleEdgesHandlingMethod;
 	};
} /* namespace NetworKit */
#endif /* FASTKONECTGRAPHREADER_H_ */
