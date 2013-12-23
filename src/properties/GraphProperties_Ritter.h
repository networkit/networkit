/*
 * GraphProperties_Ritter.h
 *
 *  Created on: 24.11.2013
 *      Author: mritter
 */

#ifndef GRAPHPROPERTIES_RITTER_H_
#define GRAPHPROPERTIES_RITTER_H_


#include "../graph/Graph.h"
#include "../io/METISGraphReader.h"
#include "ApproximateClusteringCoefficient_Hoske.h"

namespace NetworKit {

/**
 * Collection of methods for basic network properties.
 */
class GraphProperties_Ritter {
public:
	static const count INF_DIST;

	GraphProperties_Ritter();
	virtual ~GraphProperties_Ritter();

	static count diameter(const Graph& G);

	static count diameterOfTree(const Graph& T, node root);

	static count eccentricity(const Graph& G, node v);

	static std::pair<count, node> eccentricityWithNode(const Graph& G, node v);

	static Graph spanningTree(const Graph& G, node root);
};

} /* namespace NetworKit */
#endif /* GRAPHPROPERTIES_RITTER_H_ */
