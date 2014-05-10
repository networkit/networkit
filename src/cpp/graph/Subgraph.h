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
 * Methods for creating subgraphs.
 */
class Subgraph {
public:
	Subgraph();

	virtual ~Subgraph();

	/**
	 * Create a subgraph induced by a set of nodes.
	 * The returned graph G' is isomorphic (structurally identical) to the subgraph in G,
	 * but node indices are not preserved.
	 */
	static Graph fromNodes(const Graph& G, const std::unordered_set<node>& nodes);

	/** only to be used by cython - this eliminates an unnecessary copy */
	Graph* _fromNodes(const Graph& G, const std::unordered_set<node>& nodes) {
		return new Graph{std::move(fromNodes(G, nodes))};
	};
};


} /* namespace NetworKit */
#endif /* SUBGRAPH_H_ */
