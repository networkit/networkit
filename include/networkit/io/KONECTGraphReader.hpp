/*
 * KONECTGraphReader.h
 *
 *  Created on: 11.05.2018
 *      Author: Roman Bange
 */

#ifndef KONECTGRAPHREADER_H_
#define KONECTGRAPHREADER_H_

#include <unordered_map>

#include "../graph/Graph.hpp"
#include "GraphReader.hpp"

namespace NetworKit {
  class KONECTGraphReader : public GraphReader{

	public:

		/*
		* If the input graph has multiple edges, you can specify on how these edges are handled.
		* Keep in mind that NetworKit node id's start with 0 while most KONECT graphs start with 1.
		* See GraphReader.h for a closer description of the paramters.
		*
		* @param[in]	remapNodes	specifies whether node ids should be remapped if non consecutive
		* @param[in]	handlingmethod	specifies how multiple edges should be handled (only relevant if graph with multiple edges is given)
		*/
		KONECTGraphReader(bool remapNodes = false, MultipleEdgesHandling handlingmethod = DISCARD_EDGES);

		/**
		 * Given the path of an input file, read the graph contained.
		 *
		 * @param[in]	path	input file path
		 */
		virtual Graph read(const std::string& path) override;

	protected:

		bool remapNodes;
		MultipleEdgesHandling multipleEdgesHandlingMethod;
 	};
} /* namespace NetworKit */
#endif /* KONECTGRAPHREADER_H_ */
