/*
 * Graph2.h
 *
 *  Created on: 04.02.2013
 *      Author: cls
 */

#ifndef GRAPH2_H_
#define GRAPH2_H_

#include <cassert>
#include <vector>
#include <cinttypes>
#include <string>

#include "../aux/Log.h"

#define none -1

namespace EnsembleClustering {

/** Typedefs **/

typedef int64_t index; // more expressive name for an index into an array
typedef int64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based
typedef float edgeweight; // edge weight type
template<typename T> using nodemap = std::vector<T>; // more expressive name for container that is indexed by a node
template<typename T> using edgemap = std::vector<std::vector<T> >; // more expressive name for an edge data structure

class Graph2 {

protected:

	// defaults
	edgeweight defaultEdgeWeight = 1.0;
	edgeweight nullWeight = 0.0;

	// scalars
	count n; //!< number of nodes
	node maxn; //!< maximum node id / upper bound of node range

	// per node data
	nodemap<count> deg; //!< degree of each node

	// per edge data
	std::vector<std::vector<node> > adja; //!< neighbors/adjacencies
	std::vector<std::vector<edgeweight> > eweights; //!< edge weights

	/**
	 * Return the index of v in the adjacency array of u.
	 */
	index find(node u, node v) const;

public:

	Graph2(count n);

	virtual ~Graph2();

	/**
	 * Insert an undirected edge between two nodes.
	 */
	void insertEdge(node u, node v);

	/**
	 * Check if undirected edge {u,v} exists in G
	 *
	 */
	bool hasEdge(node u, node v) const;

	/**
	 * Remove undirected edge between two nodes.
	 */
	void removeEdge(node u, node v);

	/**
	 * @return Number of neighbors.
	 */
	count degree(node v) const;

	/**
	 * @return Weighted degree of @a v.
	 */
	edgeweight weightedDegree(node v) const;

	/**
	 * Return edge weight.
	 *
	 * Return 0 if edge does not exist.
	 */
	edgeweight weight(node u, node v) const;

	/**
	 * Set the weight of an edge
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void setWeight(node u, node v, edgeweight w);

	/**
	 * Add a new node to the graph and return it.
	 */
	node addNode();

	/**
	 * After calling this, nodes 1..n exist in the graph.
	 */
	void extendNodeRange(int64_t n);

	/**
	 * Return true if graph contains no nodes.
	 */
	bool isEmpty();

	/**
	 * Return the number of nodes in the graph.
	 *
	 */
	count numberOfNodes() const;

	/**
	 * Return the number of edges in the graph.
	 *
	 * This involves calculation, so store result in a if needed multiple times.
	 */
	count numberOfEdges() const;

	// TODO: lambda iterators

	/**
	 * Iterate over all nodes of the graph and execute handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle);

	/**
	 * Iterate over all nodes of the graph and execute handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle) const;

	/**
	 * Iterate over all undirected pairs of nodesand execute handler (lambda closure).
	 */
	template<typename L> void forNodePairs(L handle);

	/**
	 * Iterate in parallel over all nodes of the graph and execute handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes of the graph and execute handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle) const;

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> float parallelSumForNodes(L handle, float sum);

	/**
	 * Iterate over all edges of the graph and execute handler (lambda closure).
	 */
	template<typename L> void forEdges(L handle);

	/**
	 * Iterate in parallel over all edges of the graph and execute handler (lambda closure).
	 */
	template<typename L> void parallelForEdges(L handle);

	/**
	 * Iterate over all neighbors of a node and execute callback handletion (lamdba closure).
	 */
	template<typename L> void forallNeighborsOf(node v, L handle);

	// TODO: const iterators
};

} /* namespace EnsembleClustering */

template<typename L>
inline void EnsembleClustering::Graph2::forallNeighborsOf(node v, L handle) {
	for (node u : this->adja[v]) {
		handle(v, u);
	}
}

template<typename L>
inline void EnsembleClustering::Graph2::forNodes(L handle) {
	// TODO: this assumes that no nodes are deleted
	for (node v = 0; v < n; ++v) {
		handle(v);
	}
}

template<typename L>
inline void EnsembleClustering::Graph2::forNodes(L handle) const {
	// TODO: this assumes that no nodes are deleted
	for (node v = 0; v < n; ++v) {
		handle(v);
	}
}

template<typename L>
inline void EnsembleClustering::Graph2::parallelForNodes(L handle) {
#pragma omp parallel for
	for (node v = 0; v < n; ++v) {
		// call here
		handle(v);
	}
}

template<typename L>
inline void EnsembleClustering::Graph2::parallelForNodes(L handle) const {
#pragma omp parallel for
	for (node v = 0; v < n; ++v) {
		// call here
		handle(v);
	}
}

template<typename L>
inline float EnsembleClustering::Graph2::parallelSumForNodes(L handle, float sum) {
	sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < n; ++v) {
		// call here
		sum += handle(v);
	}
	return sum;
}

template<typename L>
inline void EnsembleClustering::Graph2::forEdges(L handle) {
	for (node u = 0; u < n; ++u) {
		for (node v : this->adja[u]) {
			if (u < v) { // {u, v} instead of (u, v); if v == -1, u < v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph2::parallelForEdges(L handle) {
#pragma omp parallel for
	for (node u = 0; u < n; ++u) {
		for (node v : this->adja[u]) {
			if (u < v) { // {u, v} instead of (u, v); if v == -1, u < v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph2::forNodePairs(L handle) {
	for (node u = 0; u < n; ++u) {
		for (node v = u + 1; v < n; ++v) {
			// call node pair function
			handle(u, v);
		}
	}
}

#endif /* GRAPH2_H_ */
