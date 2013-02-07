/*
 * Graph2.h
 *
 *  Created on: 04.02.2013
 *      Author: cls
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <cassert>
#include <vector>
#include <cinttypes>
#include <string>
#include <queue>
#include <stdexcept>

#include "../aux/Log.h"

#define none -1

namespace EnsembleClustering {

/** Typedefs **/

typedef int64_t index; // more expressive name for an index into an array
typedef int64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based
typedef float edgeweight; // edge weight type
template<typename T> using nodemap = std::vector<T>; // more expressive name for container that is indexed by a node
template<typename T> using edgemap = std::vector<std::vector<T> >;// more expressive name for an edge data structure

class Graph {

protected:

	// defaults
	static constexpr double defaultEdgeWeight = 1.0;static constexpr edgeweight nullWeight = 0.0;

	// scalars
	count n; //!< number of nodes

	// per node data
	std::vector<count> deg; //!< degree of each node

	// per edge data
	std::vector<std::vector<node> > adja; //!< neighbors/adjacencies
	std::vector<std::vector<edgeweight> > eweights; //!< edge weights

	// graph attributes
	std::string name;

	/**
	 * Return the index of v in the adjacency array of u.
	 */
	index find(node u, node v) const;

public:

	Graph(count n);

	virtual ~Graph();

	/**
	 * Set name of graph.
	 */
	void setName(std::string name);

	/*
	 * @return name of graph
	 */
	std::string getName();

	/**
	 * Get string representation
	 */
	std::string toString();

	/**
	 * Insert an undirected edge between two nodes.
	 */
	void insertEdge(node u, node v, edgeweight weight = defaultEdgeWeight);

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
	 * DEPRECATED
	 *
	 * Set self-loop weight
	 */
	void setWeight(node u, edgeweight w);

	/**
	 * @return sum of all edge weights
	 */
	edgeweight totalEdgeWeight();

	/**
	 * DEPRECATED - TODO: update clients
	 */
	edgeweight totalNodeWeight();

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
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle);

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle) const;

	/**
	 * Iterate over all undirected pairs of nodesand call handler (lambda closure).
	 */
	template<typename L> void forNodePairs(L handle);

	/**
	 * Iterate over all undirected pairs of nodesand call handler (lambda closure).
	 */
	template<typename L> void forNodePairs(L handle) const;

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodePairs(L handle);

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodePairs(L handle) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle) const;

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> float parallelSumForNodes(L handle, float sum);

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> float parallelSumForNodes(L handle, float sum) const;

	/**
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the handler
	 */
	template<typename L> float parallelSumForEdgesAndTheirWeights(L handle, float sum) const;

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forEdges(L handle);

	/**
	 * Iterate over all edges and their weights of the graph and call handler (lambda closure).
	 */
	template<typename L> void forEdgesAndTheirWeights(L handle) const;

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forEdges(L handle) const;

	/**
	 * Iterate in parallel over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForEdges(L handle);

	/**
	 * Iterate in parallel over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForEdges(L handle) const;

	/**
	 * Iterate over all neighbors of a node and call handler (lamdba closure).
	 */
	template<typename L> void forNeighborsOf(node u, L handle);

	/**
	 * Iterate over all neighbors of a node and call handler (lamdba closure).
	 */
	template<typename L> void forNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all edge weights of a node and call handler (lamdba closure).
	 */
	template<typename L> void forNeighborWeightsOf(node u, L handle) const;

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 */
	template<typename L> void forEdgesOf(node v, L handle);

	/**
	 * Iterate over nodes in breadth-first search order starting from r until connected component
	 * of r has been visited.
	 */
	template<typename L> void breadthFirstNodesFrom(node r, L handle);

	/**
	 * Iterate over edges in breadth-first search order starting from node r until connected component
	 * of r has been visited.
	 */
	template<typename L> void breadthFirstEdgesFrom(node r, L handle);

	// TODO: const iterators

	/** Adapters **/

	/**
	 * DEPRECATED - only used until the old STINGER-based interface is replaced
	 * @return the self-loop weight of this node
	 *
	 */
	edgeweight weight(node v);

};

} /* namespace EnsembleClustering */

template<typename L>
inline void EnsembleClustering::Graph::forNeighborsOf(node u, L handle) {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(v);
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::forNeighborsOf(node u, L handle) const {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(v);
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::forNeighborWeightsOf(node u, L handle) const {
	for (index i = 0; i < adja[u].size(); ++i) {
		node v = adja[u][i];
		if (v != none) {
			edgeweight ew = eweights[u][i];
			handle(v, ew);
			assert(ew == weight(u, v));
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::forNodes(L handle) {
	// TODO: this assumes that no nodes are deleted
	for (node v = 0; v < n; ++v) {
		handle(v);
	}
}

template<typename L>
inline void EnsembleClustering::Graph::forNodes(L handle) const {
	// TODO: this assumes that no nodes are deleted
	for (node v = 0; v < n; ++v) {
		handle(v);
	}
}

template<typename L>
inline void EnsembleClustering::Graph::parallelForNodes(L handle) {
#pragma omp parallel for
	for (node v = 0; v < n; ++v) {
		// call here
		handle(v);
	}
}

template<typename L>
inline void EnsembleClustering::Graph::parallelForNodes(L handle) const {
#pragma omp parallel for
	for (node v = 0; v < n; ++v) {
		// call here
		handle(v);
	}
}

template<typename L>
inline float EnsembleClustering::Graph::parallelSumForNodes(L handle,
		float sum) {
	sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < n; ++v) {
		// call here
		sum += handle(v);
	}
	return sum;
}

template<typename L>
inline float EnsembleClustering::Graph::parallelSumForNodes(L handle,
		float sum) const {
	sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < n; ++v) {
		// call here
		sum += handle(v);
	}
	return sum;
}

template<typename L>
float EnsembleClustering::Graph::parallelSumForEdgesAndTheirWeights(L handle, float sum) const {
	sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (node u = 0; u < n; ++u) {
		for (index i = 0; i < this->adja[u].size(); ++i) {
			node v = this->adja[u][i];
			edgeweight ew = this->eweights[u][i];
			if (u <= v) { // {u, v} instead of (u, v); if v == -1, u < v is not fulfilled
				sum += handle(u, v, ew);
			}
		}
	}
	return sum;
}


template<typename L>
inline void EnsembleClustering::Graph::forEdges(L handle) {
	for (node u = 0; u < n; ++u) {
		for (node v : this->adja[u]) {
			if (u <= v) { // {u, v} instead of (u, v); if v == -1, u < v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::forEdges(L handle) const {
	for (node u = 0; u < n; ++u) {
		for (node v : this->adja[u]) {
			if (u <= v) { // {u, v} instead of (u, v); if v == -1, u < v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
void EnsembleClustering::Graph::forEdgesAndTheirWeights(L handle) const {
	for (node u = 0; u < n; ++u) {
		for (index i = 0; i < this->adja[u].size(); ++i) {
			node v = this->adja[u][i];
			if (u <= v) { // {u, v} instead of (u, v); if v == -1, u < v is not fulfilled
				edgeweight ew = this->eweights[u][i];
				handle(u, v, ew);
			}
		}
	}
}


template<typename L>
inline void EnsembleClustering::Graph::parallelForEdges(L handle) {
#pragma omp parallel for
	for (node u = 0; u < n; ++u) {
		for (node v : this->adja[u]) {
			if (u <= v) { // {u, v} instead of (u, v); if v == -1, u < v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::parallelForEdges(L handle) const {
#pragma omp parallel for
	for (node u = 0; u < n; ++u) {
		for (node v : this->adja[u]) {
			if (u <= v) { // {u, v} instead of (u, v); if v == -1, u < v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::forNodePairs(L handle) {
	for (node u = 0; u < n; ++u) {
		for (node v = u + 1; v < n; ++v) {
			// call node pair function
			handle(u, v);
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::forNodePairs(L handle) const {
	for (node u = 0; u < n; ++u) {
		for (node v = u + 1; v < n; ++v) {
			// call node pair function
			handle(u, v);
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::breadthFirstNodesFrom(node r, L handle) {
	std::queue<node> q;
	count n = this->numberOfNodes();
	bool marked[n];
	for (index i = 0; i < n; ++i) {
		marked[i] = false;
	}
	q.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = q.front();
		q.pop();
		// apply function
		handle(u);
		this->forNeighborsOf(u, [&](node v) {
			// filtering edges is not necessary because only out-edges are considered by stinger
				if (!marked[v]) {
					TRACE("pushing node " << v);
					q.push(v);
					marked[v] = true;
				}
			});
	} while (!q.empty());
}

template<typename L>
inline void EnsembleClustering::Graph::forEdgesOf(node u, L handle) {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(u, v);
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::parallelForNodePairs(L handle) {
#pragma omp parallel for
	for (node u = 0; u < n; ++u) {
		for (node v = u + 1; v < n; ++v) {
			// call node pair function
			handle(u, v);
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::parallelForNodePairs(L handle) const {
#pragma omp parallel for
	for (node u = 0; u < n; ++u) {
		for (node v = u + 1; v < n; ++v) {
			// call node pair function
			handle(u, v);
		}
	}
}

template<typename L>
inline void EnsembleClustering::Graph::breadthFirstEdgesFrom(node r, L handle) {
	// TODO:
	throw std::runtime_error("TODO");
}

#endif /* GRAPH_H_ */
