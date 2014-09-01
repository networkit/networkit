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
 * @ingroup community
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
	 *
	 * @param G The graph.
	 * @param zeta A community clustering of @a G.
	 */
	virtual void run(const Graph& G, const Partition& zeta);

	/**
	 * Returns the coarsened Graph.
	 * @return The coarsened graph.
	 */
	virtual Graph getGraph();

	/** 
	 * Maps community id to node id in the community graph.
	 *
	 * @return Map containing community id to node id mappings.
	 */
	virtual std::map<index, node> getCommunityToNodeMap();

	/** 
	 * Maps node id in the community graph to community id.
	 *
	 * @return Map containing node id to community id mappins.
	 */
	virtual std::map<node, index> getNodeToCommunityMap(); 


private:
	Graph Gcom;		//!< the current community graph
	std::vector<node> communityToSuperNode;		//!< community id -> node id in community graph

};

} /* namespace NetworKit */
#endif /* COMMUNITYGRAPH_H_ */
