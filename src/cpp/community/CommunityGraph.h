/*
 * CommunityGraph.h
 *
 *  Created on: 22.08.2013
 *      Author: cls
 */

#ifndef COMMUNITYGRAPH_H_
#define COMMUNITYGRAPH_H_

#include "../structures/Partition.h"

namespace NetworKit {

/**
 * The CommunityGraph class represents a Graph coarsened according to communities. Each node in the CommunityGraph
 * represents a community. Edge weights are the weights of inter-community cuts.
 */
class CommunityGraph {
public:
	/** Default destructor */
	virtual ~CommunityGraph();

	/**
	 * Creates a coarsened graph of @a G according to communities in @a zeta. Edge weights are the weights of
	 * inter-community cuts.
	 */
	virtual void run(const Graph& G, const Partition& zeta);

	/**
	 * Returns the coarsened Graph.
	 * @return The coarsened graph.
	 */
	virtual Graph getGraph();

	/** only to be used by cython - this eliminates an unnecessary copy */
	Graph* _getGraph() {
		return new Graph{std::move(getGraph())};
	};

	/** 
	 * Maps community id to id of node in the community graph.
	 */
	virtual std::map<index, node> getCommunityToNodeMap();

	/** 
	 * Maps node id in the community graph to community id.
	 */
	virtual std::map<node, index> getNodeToCommunityMap(); 


private:
	Graph Gcom;		//!< the current community graph
	std::vector<node> communityToSuperNode;		//!< community id -> node id in community graph

};

} /* namespace NetworKit */
#endif /* COMMUNITYGRAPH_H_ */
