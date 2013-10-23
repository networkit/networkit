/*
 * CommunityGraph.h
 *
 *  Created on: 22.08.2013
 *      Author: cls
 */

#ifndef COMMUNITYGRAPH_H_
#define COMMUNITYGRAPH_H_

#include "../graph/Graph.h"
#include "../clustering/Clustering.h"

namespace NetworKit {

class CommunityGraph {
public:

	CommunityGraph();

	virtual ~CommunityGraph();

	/**
	 * Create a graph coarsened according to communities. Edge weights are the weights of
	 * inter-community cuts.
	 */
	virtual Graph get(const Graph& G, const Clustering& zeta);
};

} /* namespace NetworKit */
#endif /* COMMUNITYGRAPH_H_ */
