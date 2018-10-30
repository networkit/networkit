/*
 * KONECTGraphReader.h
 *
 *  Created on: 11.05.2018
 *      Author: Roman Bange
 */

#ifndef KONECTGRAPHREADER_H_
#define KONECTGRAPHREADER_H_

#include <unordered_map>

#include "../graph/Graph.h"
#include "GraphReader.h"

namespace NetworKit {
  class KONECTGraphReader : public GraphReader{

	public:
		[[deprecated("New Constructor - Use standardized seperator and Graph.removeSelfLoops instead")]]
		KONECTGraphReader(char separator, bool ignoreLoops=false) : KONECTGraphReader(false, DISCARD_EDGES){
			if(separator != ' ' && separator != '\t')
				throw std::runtime_error("Separator is not supported anymore. Use space or tab as standardized for KONECT formats");
			if (ignoreLoops)
				ERROR("IgnoreLoops not supported anymore - use Graph.removeSelfLoops after reading.");
		};

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
