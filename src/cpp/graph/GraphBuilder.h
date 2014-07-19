/*
 * GraphBuilder.h
 *
 *  Created on: 15.07.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef GRAPH_BUILDER_H
#define GRAPH_BUILDER_H

#include <vector>

#include "../Globals.h"
#include "Graph.h"

namespace NetworKit {

class GraphBuilder {
protected:
	count n; //!< current number of nodes

	bool weighted; //!< true if the graph will be weighted, false otherwise
	bool directed; //!< true if the graph will be directed, false otherwise

	std::vector< std::vector<node> > halfEdges;
	std::vector< std::vector<edgeweight> > halfEdgeWeights;

	index indexHalfEdgeArray(node u, node v) const;

public:
	GraphBuilder(count n = 0, bool weighted = false, bool directed = false);

	bool isWeighted() const { return weighted; }
	bool isDirected() const { return directed; }
	count numberOfNodes() const { return n; }

	node addNode();
	void addEdge(node u, node v, edgeweight ew = defaultEdgeWeight);
	void increaseWeight(node u, node v, edgeweight ew);

	Graph toGraph(bool parallel = true) { return parallel ? toGraphParallel() : toGraphSequential(); }

	/**
	 * Iterate over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void forNodes(L handle) const;

	/**
	 * Iterate randomly over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void parallelForNodes(L handle) const;

	/**
	 * Iterate over all undirected pairs of nodes and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void forNodePairs(L handle) const;


	/**
	 * Iterate over all undirected pairs of nodes in parallel and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void parallelForNodePairs(L handle) const;

private:
	Graph toGraphParallel();
	Graph toGraphSequential();
};

template<typename L>
void GraphBuilder::forNodes(L handle) const {
	for (node v = 0; v < n; v++) {
		handle(v);
	}
}

template<typename L>
void GraphBuilder::parallelForNodes(L handle) const {
	#pragma omp parallel for schedule(dynamic)
	for (node v = 0; v < n; v++) {
		handle(v);
	}
}

template<typename L>
void GraphBuilder::forNodePairs(L handle) const {
	for (node u = 0; u < n; u++) {
		for (node v = u + 1; v < n; v++) {
			handle(u, v);
		}
	}
}

template<typename L>
void GraphBuilder::parallelForNodePairs(L handle) const {
	#pragma omp parallel for schedule(dynamic)
	for (node u = 0; u < n; u++) {
		for (node v = u + 1; v < n; v++) {
			handle(u, v);
		}
	}
}

} /* namespace NetworKit */

#endif /* GRAPH_BUILDER_H */
