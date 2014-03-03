/*
 * Graph.h
 *
 *  Created on: 04.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <functional>
#include <cassert>
#include <vector>
#include <cinttypes>
#include <string>
#include <queue>
#include <stdexcept>
#include <map>
#include <set>
#include <sstream>
#include <limits>
#include <cstdint>
#include <algorithm>
// #include <tbb/concurrent_vector.h>

#include "Coordinates.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Debug.h"
#include "../auxiliary/PriorityQueue.h"
#include "../Globals.h"
#include "../viz/Point.h"


namespace NetworKit {

/** Typedefs **/

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based
typedef double edgeweight; // edge weight type

constexpr index none = std::numeric_limits<index>::max();

//#define Vector std::vector // TODO: test tbb::concurrent_vector
template <typename T> using Vector = std::vector<T>;

/**
 * An undirected graph (with optional weights) and parallel iterator methods.
 */
class Graph {

protected:

	// scalars
	count n; //!< current number of nodes
	count m; //!< current number of edges
	node z; //!< current upper bound of node ids
	count t; //!< current time step
	bool weighted; //!< true if this graph has been marked as weighted.

	// per node data
	Vector<count> deg; //!< degree of each node (size of neighborhood)
	Vector<bool> exists; //!< exists[v] is true if node v has not been removed from the graph
	Coordinates<float> coordinates; //!< coordinates of nodes (if present)

	// per edge data
	Vector<Vector<node> > adja; //!< neighbors/adjacencies
	Vector<Vector<edgeweight> > eweights; //!< edge weights

	// graph attributes
	std::string name;

	// user-defined edge attributes

	//	attribute maps storage

	std::vector<std::vector<std::vector<double> > > edgeMaps_double; // contains edge maps (u, v) -> double

	// defaults

	std::vector<double> edgeAttrDefaults_double; // stores default value for edgeMaps_double[i] at index i


	/**
	 * Return the index of v in the adjacency array of u.
	 */
	index find(node u, node v) const;

public:

	// defaults
	static constexpr double defaultEdgeWeight = 1.00;
	static constexpr edgeweight nullWeight = 0.0;

	/** ATTRIBUTE ABSTRACT BASE CLASSES **/

	class NodeAttribute {
		// abstract
	};

	class EdgeAttribute {
		// abstract
	};


	/** GRAPH INTERFACE **/

	Graph();

	/** 
	 * Create a graph of n nodes.
	 */
	Graph(count n);

	Graph(const Graph& other) = default;

	Graph(Graph&& other) = default;

	virtual ~Graph();

	Graph& operator=(Graph&& other) = default;

	Graph& operator=(const Graph& other) = default;



	/**
	 * Set name of graph.
	 */
	void setName(std::string name);

	/*
	 * @return name of graph
	 */
	std::string getName();

	/**
	 * Mark this graph as a weighted graph, i.e. as a graph containing edges
	 * with weight other than 1.0.
	 * The insertion of edges with weights other than 1.0 does not automatically
	 * mark the graph as weighted.
	 */
	void markAsWeighted();


	/**
	 * Return if this graph has been marked as a weighted graph.
	 */
	bool isMarkedAsWeighted() const;

	/**
	 * Get string representation
	 */
	std::string toString();

	/**
	 * Insert an undirected edge between two nodes.
	 */
	void addEdge(node u, node v, edgeweight weight = defaultEdgeWeight);

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
	 * Merges edge {u,v} to become a supernode. Edges to u and v are
	 * rewired, multiple edges merged and their weights added.
	 * The vertex weights of @a u and @a v are added.
	 * A self-loop is only created if @a discardSelfLoop is set to false.
	 *
	 * @return New node that has been created if u != v. Otherwise none.
	 */
	node mergeEdge(node u, node v, bool discardSelfLoop = true);

	/**
	 * @return Number of neighbors.
	 */
	count degree(node v) const;

	/**
	 * @return Smallest neighborhood size (does not have to be unique).
	 */
	count minDegree() const;

	/**
	 * @return Index of vertex with smallest neighborhood size (does not have to be
	 * unique).
	 */
	index argminDegree() const;

	/**
	 * @return Largest neighborhood size (does not have to be unique).
	 */
	count maxDegree() const;

	/**
	 * @return Index of vertex with largest neighborhood size (does not have to be
	 * unique).
	 */
	index argmaxDegree() const;

	/**
	 * @return Weighted degree of @a v.
	 */
	edgeweight weightedDegree(node v) const;

	/**
	 * @return Volume of the node, which is the
	 * weighted degree with self-loops counted twice.
	 */
	edgeweight volume(node v) const;

	/**
	 * @return Random (uuid) neighbor of @a v. None if degree is zero.
	 */
	node randomNeighbor(node v) const;



	/** EDGE ATTRIBUTE GETTERS **/

	/**
	 * Return edge weight.
	 *
	 * Return 0 if edge does not exist.
	 */
	edgeweight weight(node u, node v) const;

	/**
	 * @return attribute of type double for an edge.
	 *
	 * @param[in]	u	node
	 * @param[in]	v	node
	 * @param[in]	attrId	attribute id
	 */
	double attribute_double(node u, node v, int attrId) const;

	/**  EDGE ATTRIBUTE SETTERS */

	/**
	 * Set the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void setWeight(node u, node v, edgeweight w);


	/**
	 * Increase the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void increaseWeight(node u, node v, edgeweight w);

	/**
	 * Set edge attribute of type double If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	attr	double edge attribute
	 */
	void setAttribute_double(node u, node v, int attrId, double attr);

	/** SUMS **/

	/**
	 * @return sum of all edge weights
	 */
	edgeweight totalEdgeWeight() const;


	/** NODE MODIFIERS **/

	/**
	 * Add a new node to the graph and return it.
	 */
	node addNode();

	/**
	 * Add a new node to the graph with coordinates @a x and @y and return it.
	 */
	node addNode(float x, float y);

	/**
	 * Remove an isolated node from the graph.
	 *
	 * Although it would be convenient to remove all incident edges at the same time,
	 * this causes complications for dynamic applications. Therefore, removeNode is an
	 * atomic event. All incident edges need to be removed first and an exception is thrown
	 * otherwise.
	 */
	void removeNode(node u);

	/**
	 * Check if node exists in the graph.
	 */
	bool hasNode(node u) const;


	/** GLOBAL PROPERTIES **/

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
	 *	 */
	count numberOfEdges() const;

	/**
	 * @return the number of loops {v, v} in the graph.
	 *
	 * This involves calculation, so store result if needed multiple times.
	 */
	count numberOfSelfLoops() const;


	/**
	 * Get an upper bound for the node ids in the graph.
	 */
	index upperNodeIdBound() const;

	/** DYNAMICS **/

	/**
	 * Trigger a time step - increments counter.
	 */
	void timeStep();

	/**
	 * Get time step counter.
	 */
	count time();


	/** COORDINATES **/

	void setCoordinate(node v, Point<float> value) {
		coordinates.setCoordinate(v, value);
	}

	Point<float>& getCoordinate(node v) {
		return coordinates.getCoordinate(v);
	}

	float minCoordinate(count dim) {
		return coordinates.minCoordinate(dim);
	}

	float maxCoordinate(count dim) {
		return coordinates.maxCoordinate(dim);
	}

	void initCoordinates() {
		coordinates.init(z);
	}

	/** ATTRIBUTES **/

	/**
	 * Add new edge map for an attribute of type double.
	 */
	int addEdgeAttribute_double(double defaultValue);

	/** NODE ITERATORS **/

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle);

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle) const;

	/**
	 * Iterate randomly over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodesInRandomOrder(L handle);

	/**
	 * Iterate randomly over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodesInRandomOrder(L handle) const;

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	template<typename C, typename L> void forNodesWhile(C condition, L handle);

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	template<typename C, typename L> void forNodes(C condition, L handle) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 */
	template<typename L> void balancedParallelForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 */
	template<typename L> void balancedParallelForNodes(L handle) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle) const;

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
	 * Iterate over nodes in breadth-first search order starting from r until connected component
	 * of r has been visited.
	 */
	template<typename L> void breadthFirstNodesFrom(node r,
			std::vector<int>& marked, L handle);

	/**
	 * Iterate over edges in breadth-first search order starting from node r until connected component
	 * of r has been visited.
	 */
	template<typename L> void breadthFirstEdgesFrom(node r, L handle);

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 *
	 * @param[in]	attrKey		attribute key
	 * @param[in]	handle		takes parameters (v, a) where a is a node attribute
	 */
	template<typename L> void forNodesWithAttribute(std::string attrKey,
			L handle);

	/** EDGE ITERATORS **/

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forEdges(L handle);

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
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 *
	 */
	template<typename L> void forWeightedEdges(L handle);

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 */
	template<typename L> void forWeightedEdges(L handle) const;

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 *
	 */
	template<typename L> void parallelForWeightedEdges(L handle);

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 */
	template<typename L> void parallelForWeightedEdges(L handle) const;

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 *	@param[in]	attrId		attribute id
	 *	@param[in]	handle 		takes arguments (u, v, a) where a is an edge attribute of edge {u, v}
	 *
	 */
	template<typename L> void forEdgesWithAttribute_double(int attrId,
			L handle);

	/**
	 * Iterate over all edges of the const graph and call handler (lambda closure).
	 *
	 *	@param[in]	attrId		attribute id
	 *	@param[in]	handle 		takes arguments (u, v, a) where a is an edge attribute of edge {u, v}
	 *
	 */
	template<typename L> void forEdgesWithAttribute_double(int attrId,
			L handle) const;

	/** NEIGHBORHOOD ITERATORS **/

	/**
	 * Iterate over all neighbors of a node and call handler (lamdba closure).
	 *
	 * (Note that a node is its own neighbor if there is a self-loop.)
	 */
	template<typename L> void forNeighborsOf(node u, L handle);

	/**
	 * Iterate over all neighbors of a node and call handler (lamdba closure).
	 */
	template<typename L> void forNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all edge weights of a node and call handler (lamdba closure).
	 */
	template<typename L> void forWeightedNeighborsOf(node u, L handle);

	/**
	 * Iterate over all edge weights of a node and call handler (lamdba closure).
	 */
	template<typename L> void forWeightedNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 */
	template<typename L> void forEdgesOf(node u, L handle);

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 */
	template<typename L> void forEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all incident edges of a node in neighborhood-size-increasing order and call handler (lamdba closure).
	 */
	template<typename L> void forEdgesOfInDegreeIncreasingOrder(node u, L handle);

	/**
	 * Iterate over all incident edges of a node in neighborhood-size-increasing order and call handler (lamdba closure).
	 */
	template<typename L> void forEdgesOfInDegreeIncreasingOrder(node u, L handle) const;

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 *
	 */
	template<typename L> void forWeightedEdgesOf(node u, L handle);

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 *
	 */
	template<typename L> void forWeightedEdgesOf(node u, L handle) const;

	/** REDUCTION ITERATORS **/

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle) const;

	/**
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForWeightedEdges(L handle) const;


	/** Collections **/

	/**
	 * Return list of nodes
	 */
	std::vector<node> nodes();

	/**
	 * Return list of edges as node pairs.
	 */
	std::vector<std::pair<node, node> > edges();


	/**
	 * Return list of neighbors for given node.
	 */
	std::vector<node> neighbors(node u);


};

} /* namespace NetworKit */

template<typename L>
inline void NetworKit::Graph::forNeighborsOf(node u, L handle) {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forNeighborsOf(node u, L handle) const {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(v);
		}
	}
}


template<typename L>
inline void NetworKit::Graph::forWeightedNeighborsOf(node u, L handle) {
	for (index i = 0; i < (index) adja[u].size(); ++i) {
		node v = adja[u][i];
		if (v != none) {
			edgeweight ew = eweights[u][i];
			handle(v, ew);
			assert(ew == weight(u, v));
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forWeightedNeighborsOf(node u, L handle) const {
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
inline void NetworKit::Graph::forNodes(L handle) {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forNodes(L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForNodes(L handle) {
#pragma omp parallel for
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForNodes(L handle) const {
#pragma omp parallel for
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::balancedParallelForNodes(L handle) {
#pragma omp parallel for schedule(guided) // TODO: define min block size (and test it!)
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::balancedParallelForNodes(L handle) const {
#pragma omp parallel for schedule(guided) // TODO: define min block size (and test it!)
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline double NetworKit::Graph::parallelSumForNodes(L handle) {
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			sum += handle(v);
		}
	}
	return sum;
}

template<typename L>
inline double NetworKit::Graph::parallelSumForNodes(L handle) const {
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			sum += handle(v);
		}
	}
	return sum;
}

template<typename L>
double NetworKit::Graph::parallelSumForWeightedEdges(L handle) const {
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->adja[u].size(); ++i) {
			node v = this->adja[u][i];
			edgeweight ew = this->eweights[u][i];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				sum += handle(u, v, ew);
			}
		}
	}
	return sum;
}

template<typename L>
inline void NetworKit::Graph::forEdges(L handle) {
	for (node u = 0; u < z; ++u) {
		for (node v : this->adja[u]) {
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forEdges(L handle) const {
	for (node u = 0; u < z; ++u) {
		for (node v : this->adja[u]) {
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForEdges(L handle) {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (node v : this->adja[u]) {
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForEdges(L handle) const {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (node v : this->adja[u]) {
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forNodePairs(L handle) {
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				// call node pair function
				if (exists[v]) {
					handle(u, v);
				}
			}
		}

	}
}

template<typename L>
inline void NetworKit::Graph::forNodePairs(L handle) const {
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				// call node pair function
				if (exists[v]) {
					handle(u, v);
				}
			}
		}

	}
}

template<typename L>
inline void NetworKit::Graph::breadthFirstNodesFrom(node r,
		std::vector<int>& marked, L handle) {
	std::queue<node> q;
	q.push(r); // enqueue root
	marked[r] = 1;
	do {
		node u = q.front();
		q.pop();
		// apply function
		handle(u);
		this->forNeighborsOf(u, [&](node v) {
			if (marked[v] == 0) {
				q.push(v);
				marked[v] = 1;
			}
		});
	} while (!q.empty());
}

template<typename L>
inline void NetworKit::Graph::forEdgesOf(node u, L handle) {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(u, v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forEdgesOf(node u, L handle) const {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(u, v);
		}
	}
}

template<typename L>
void NetworKit::Graph::forEdgesOfInDegreeIncreasingOrder(node u, L handle) const {
	// TODO: iterating over neighbors ordered by degree does not need privileged access to graphs data structure and should
	// therefore be implemented inside the algorithm. - cls
	auto hasSmallerDegree = [&](node v1, node v2) {
		return degree(v1) < degree(v2); // FIXME
	};

	Vector<node> neighbors = adja[u];
	std::sort(neighbors.begin(), neighbors.end(), hasSmallerDegree);

	for (node v : neighbors) {
		if (v != none) {
			handle(u, v);
		}
	}
}

template<typename L>
void NetworKit::Graph::forEdgesOfInDegreeIncreasingOrder(node u, L handle) {
	auto hasSmallerDegree = [&](node v1, node v2) {
		return degree(v1) < degree(v2); // FIXME
	};

	Vector<node> neighbors = adja[u];
	std::sort(neighbors.begin(), neighbors.end(), hasSmallerDegree);

	for (node v : neighbors) {
		if (v != none) {
			handle(u, v);
		}
	}
}


template<typename L>
inline void NetworKit::Graph::parallelForNodePairs(L handle) {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				// call node pair function
				if (exists[v]) {
					handle(u, v);
				}
			}
		}

	}
}

template<typename L>
inline void NetworKit::Graph::parallelForNodePairs(L handle) const {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				// call node pair function
				if (exists[v]) {
					handle(u, v);
				}
			}
		}

	}
}

template<typename L>
inline void NetworKit::Graph::breadthFirstEdgesFrom(node r, L handle) {
	// TODO: implement BFS iterator for edges
	throw std::runtime_error("TODO");
}

template<typename L>
inline void NetworKit::Graph::forWeightedEdges(L handle) {
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				edgeweight w = this->eweights[u][vi];
				handle(u, v, w);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forWeightedEdges(L handle) const {
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				edgeweight w = this->eweights[u][vi];
				handle(u, v, w);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForWeightedEdges(L handle) {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				edgeweight w = this->eweights[u][vi];
				handle(u, v, w);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForWeightedEdges(L handle) const {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				edgeweight w = this->eweights[u][vi];
				handle(u, v, w);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forWeightedEdgesOf(node u, L handle) {
	const count asize = (count) adja[u].size();
	for (index i = 0; i < asize; ++i) {
		node v = adja[u][i];
		if (v != none) {
			edgeweight ew = eweights[u][i];
			handle(u, v, ew);
//			assert(ew == weight(u, v));
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forWeightedEdgesOf(node u, L handle) const {
	for (index i = 0; i < adja[u].size(); ++i) {
		node v = adja[u][i];
		if (v != none) {
			edgeweight ew = eweights[u][i];
			handle(u, v, ew);
			assert(ew == weight(u, v));
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forNodesWithAttribute(std::string attrKey,
		L handle) {
	// get nodemap for attrKey

//	auto nodeMap; // ?
//
//	auto findIdPair = this->attrKey2IdPair.find(attrKey);
//	if (findIdPair != this->attrKey2IdPair.end()) {
//		std::pair<index, index> idPair = findIdPair->second;
//		index typeId = idPair.first;
//		index mapId = idPair.second;
//
//		// nodemaps are in a vector, one for each node attribute type int, double, NodeAttribute
//		switch (typeId) {
//		case 0:
//			nodeMap = this->nodeMapsInt[mapId];
//			break;
//		case 1:
//			nodeMap = this->nodeMapsdouble[mapId];
//			break;
//		}
//
//		// iterate over nodes and call handler with attribute
//		this->forNodes([&](node u) {
//			auto attr = nodeMap[u];
//			handle(u, attr);
//		});
//	} else {
//		throw std::runtime_error("node attribute not found");
//	}

	// TODO: iterate over nodes with atributes
	throw std::runtime_error("TODO");
}

template<typename L>
inline void NetworKit::Graph::forEdgesWithAttribute_double(int attrId,
		L handle) {
	std::vector<std::vector<double> > edgeMap = this->edgeMaps_double[attrId];
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < (index) adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			double attr = edgeMap[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v, attr);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forEdgesWithAttribute_double(int attrId,
		L handle) const {
	std::vector<std::vector<double> > edgeMap = this->edgeMaps_double[attrId];
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < (index) adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			double attr = edgeMap[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v, attr);
			}
		}
	}
}

template<typename C, typename L>
inline void NetworKit::Graph::forNodesWhile(C condition, L handle) {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (!condition()) {
				break; // if condition does not hold, break from loop and do not call handle
			}
			handle(v);

		}
	}
}



template<typename C, typename L>
inline void NetworKit::Graph::forNodes(C condition, L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (!condition()) {
				break; // if condition does not hold, break from loop and do not call handle
			}
			handle(v);
		}
	}
}

template<typename L>
void NetworKit::Graph::forNodesInRandomOrder(L handle) {
	std::vector<node> randVec(z);
	for (node v = 0; v < z; ++v) {
		randVec[v] = v;
	}
	random_shuffle(randVec.begin(), randVec.end());

	for (node v = 0; v < z; ++v) {
		node randv = randVec[v];
		if (exists[randv]) {
			handle(randv);
		}
	}
}

template<typename L>
void NetworKit::Graph::forNodesInRandomOrder(L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}


#endif /* GRAPH_H_ */
