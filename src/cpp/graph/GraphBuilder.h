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
	count selfloops; //!< currently encountered number of self loops

	bool weighted; //!< true if the graph will be weighted, false otherwise
	bool directed; //!< true if the graph will be directed, false otherwise
	bool usedirectswap; //!< true if edges are moved directly into the graph without constructing backedges

	std::vector< std::vector<node> > halfEdges;
	std::vector< std::vector<edgeweight> > halfEdgeWeights;

	index indexHalfEdgeArray(node u, node v) const;

public:
	/**
	 * Creates a new GraphBuilder. GraphBuilder supports the basic methods needed to create a new graph (addNode, addEdge, setWeight, increaseWeight). It is designed to be much faster for graph creation, but the speed comes with a restriction:
	 * For undirected graphs GraphBuilder will handle u->v and v->u as two different edges. Keep that in mind when using setWeight and increaseWeight.
	 * GraphBuilder allows parallelization in a special way. It's internal data structure saves edges only at the source node. As long as edges from node u are only added/changed by thread t1, every other thread can modifier edges not starting in u.
	 * addNode is not threadsafe.
	 * @param n Number of nodes.
	 * @param weighted If set to <code>true</code>, the graph has edge weights.
	 * @param directed If set to @c true, the graph will be directed.
	 */
	GraphBuilder(count n = 0, bool weighted = false, bool directed = false, bool directSwap = false);

	/**
	 * Returns <code>true</code> if this graph supports edge weights other than 1.0.
	 * @return <code>true</code> if this graph supports edge weights other than 1.0.
	 */
	bool isWeighted() const { return weighted; }

	/**
	 * Return <code>true</code> if this graph supports directed edges.
	 * @return </code>true</code> if this graph supports directed edges.
	 */
	bool isDirected() const { return directed; }

	/**
	 * Return <code>true</code> if this graph builder uses direct swaps.
	 * @return </code>true</code> if this graph builder uses direct swaps.
	 */
	bool useDirectSwap() const { return usedirectswap; }

	/**
	 * Return <code>true</code> if graph contains no nodes.
	 * @return <code>true</code> if graph contains no nodes.
	 */
	bool isEmpty() const { return n == 0; }

	/**
	 * Return the number of nodes in the graph.
	 * @return The number of nodes.
	 */
	count numberOfNodes() const { return n; }

 	/**
	 * Get an upper bound for the node ids in the graph.
	 * @return An upper bound for the node ids.
	 */
	index upperNodeIdBound() const { return n; }

	/**
	 * Add a new node to the graph and return it.
	 * @return The new node.
	 */
	node addNode();
	
	/**
	 * Insert an edge between the nodes @a u and @a v. If the graph is weighted you can optionally
	 * set a weight for this edge. The default weight is 1.0.
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @param weight Optional edge weight.
	 */
	void addEdge(node u, node v, edgeweight ew = defaultEdgeWeight);

	/**
	 * Set the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void setWeight(node u, node v, edgeweight ew);

	/**
	 * Increase the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void increaseWeight(node u, node v, edgeweight ew);

	/**
	 * Generates a Graph instance. The graph builder will be reseted at the end.
	 */
	Graph toGraph(bool parallel = true) {
		Graph G(n, weighted, directed);
		if (usedirectswap) directSwap(G);
		else parallel ? toGraphParallel(G) : toGraphSequential(G);
		return G;
	}

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
	void toGraphParallel(Graph &G);
	void toGraphSequential(Graph &G);
	void directSwap(Graph &G);

	void reset();
	template <typename T>
	static void copyAndClear(std::vector<T>& source, std::vector<T>& target);
	static void correctNumberOfEdges(Graph& G, count numberOfSelfLoops);
	static bool checkConsistency(Graph& G);
};

template <typename T>
void GraphBuilder::copyAndClear(std::vector<T>& source, std::vector<T>& target) {
	std::copy(source.begin(), source.end(), std::back_inserter(target));
	source.clear();	
}

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
