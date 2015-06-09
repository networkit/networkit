/*
 * Subgraph.h
 *
 *  Created on: Jun 13, 2013
 *      Author: forigem
 */

#ifndef SUBGRAPH_H_
#define SUBGRAPH_H_

#include <unordered_set>
#include <unordered_map>

#include "Graph.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Methods for creating subgraphs.
 */
class Subgraph {
public:

	/**
	 * Create a subgraph induced by the set @a nodes.
	 *
	 * @param G The graph.
	 * @param nodes A subset of nodes of @a G which induce the subgraph.
	 * @return The subgraph induced by @a nodes.
	 * @note The returned graph G' is isomorphic (structurally identical) to the subgraph in G,
	 * but node indices are not preserved.
	 */
	static Graph fromNodes(const Graph& G, const std::unordered_set<node>& nodes);
};


} /* namespace NetworKit */
#endif /* SUBGRAPH_H_ */
