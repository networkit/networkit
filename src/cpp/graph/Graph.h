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
#include <stack>
#include <stdexcept>
#include <map>
#include <set>
#include <sstream>
#include <limits>
#include <cstdint>
#include <algorithm>
// #include <tbb/concurrent_vector.h>

#include "AbstractGraph.h"
#include "Coordinates.h"
#include "../auxiliary/Log.h"
#include "../Globals.h"
#include "../viz/Point.h"

namespace NetworKit {

/**
 * An undirected graph (with optional weights) and parallel iterator methods.
 */
class Graph final : public AbstractGraph {

protected:

	// per node data
	_Vector<count> deg; //!< degree of each node (size of neighborhood)

	// per edge data
	_Vector<_Vector<node> > adja; //!< neighbors/adjacencies
	_Vector<_Vector<edgeweight> > eweights; //!< edge weights

	// user-defined edge attributes
	std::vector<std::vector<std::vector<double> > > edgeMaps_double; // contains edge maps (u, v) -> double
	std::vector<double> edgeAttrDefaults_double; // stores default value for edgeMaps_double[i] at index i

	/**
	 * Return the index of v in the adjacency array of u.
	 */
	index find(node u, node v) const;

public:

	/** ATTRIBUTE ABSTRACT BASE CLASSES **/

	// class NodeAttribute {
	// 	// abstract
	// };

	// class EdgeAttribute {
	// 	// abstract
	// };


	/** GRAPH INTERFACE **/

	/** 
	 * Create a graph of n nodes.
	 */
	Graph(count n = 0, bool weighted = false);

	Graph(const Graph& other) = default;

	Graph(Graph&& other) = default;

	virtual ~Graph();

	Graph& operator=(Graph&& other) = default;

	Graph& operator=(const Graph& other) = default;


	/**
	 * Calculate an approximation of the memory used by this graph. Only memory increasing with the
	 * number of edges or nodes of this graph is taken into account. 
	 */
	count getMemoryUsage() const ;

	/**
	 * Try to save some memory by shrinking internal data structures of the graph. Only run this
	 * once you finished editing the graph. Otherwise it will cause unnecessary reallocation of
	 * memory. 
	 */
	void shrinkToFit();

	/** Only to be used from Cython */
	void stealFrom(Graph& input);


	/** NODE MODIFIERS **/

	/**
	 * Add a new node to the graph and return it.
	 */
	node addNode() override;

	/**
	 * Add a new node to the graph with coordinates @a x and @y and return it.
	 */
	node addNode(float x, float y) override;


	/** NODE PROPERTIES **/
 
	/**
	 * @return true if the node is isolated (= degree is 0)
	 */
	bool isIsolated(node v) const;

	/**
	 * Return the number of neighbors for node v.
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
	 * @return Random neighbor of @a v. None if degree is zero.
	 */
	node randomNeighbor(node v) const;


	/** EDGE MODIFIERS **/

	/**
	 * Insert an undirected edge between two nodes.
	 */
	void addEdge(node u, node v, edgeweight weight = defaultEdgeWeight);

	/**
	 * Remove undirected edge between two nodes.
	 */
	void removeEdge(node u, node v);

	/**
	 * Check if undirected edge {u,v} exists in G
	 *
	 */
	bool hasEdge(node u, node v) const;

	/**
	 * Merges edge {u,v} to become a supernode. Edges to u and v are
	 * rewired, multiple edges merged and their weights added.
	 * The vertex weights of @a u and @a v are added.
	 * A self-loop is only created if @a discardSelfLoop is set to false.
	 *
	 * @return New node that has been created if u != v. Otherwise none.
	 */
	node mergeEdge(node u, node v, bool discardSelfLoop = true);


	/** GLOBAL PROPERTIES **/

	/** 
	 * Return true if this graph supports directed edges.
	 */
	bool isDirected() const;


	/** EDGE ATTRIBUTES **/

	/**
	 * Return edge weight.
	 *
	 * Return 0 if edge does not exist.
	 */
	edgeweight weight(node u, node v) const;

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
	 * Add new edge map for an attribute of type double.
	 */
	int addEdgeAttribute_double(double defaultValue);

	/**
	 * @return attribute of type double for an edge.
	 *
	 * @param[in]	u	node
	 * @param[in]	v	node
	 * @param[in]	attrId	attribute id
	 */
	double attribute_double(node u, node v, int attrId) const;

	/**
	 * Set edge attribute of type double If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	attr	double edge attribute
	 */
	void setAttribute_double(node u, node v, int attrId, double attr);


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
	template<typename L> void forEdgesWithAttribute_double(int attrId, L handle);

	/**
	 * Iterate over all edges of the const graph and call handler (lambda closure).
	 *
	 *	@param[in]	attrId		attribute id
	 *	@param[in]	handle 		takes arguments (u, v, a) where a is an edge attribute of edge {u, v}
	 *
	 */
	template<typename L> void forEdgesWithAttribute_double(int attrId, L handle) const;






	/** OLD STUFF **/

	/** NODE ITERATORS **/

	/**
	 * Iterate over nodes in breadth-first search order starting from r until connected component
	 * of r has been visited.
	 */
	template<typename L> void breadthFirstNodesFrom(node r, std::vector<int>& marked, L handle);

	template<typename L> void BFSfrom(node r, L handle);

	template<typename L> void DFSfrom(node r, L handle);

	/**
	 * Iterate over edges in breadth-first search order starting from node r until connected component
	 * of r has been visited.
	 */
	template<typename L> void breadthFirstEdgesFrom(node r, L handle);


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
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForWeightedEdges(L handle) const;





};

} /* namespace NetworKit */


/** EDGE ITERATORS **/

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
inline void NetworKit::Graph::forWeightedEdges(L handle) {
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (weighted) {
					edgeweight w = this->eweights[u][vi];
					handle(u, v, w);
				} else {
					handle(u, v, defaultEdgeWeight);
				}
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
				if (weighted) {
					edgeweight w = this->eweights[u][vi];
					handle(u, v, w);
				} else {
					handle(u, v, defaultEdgeWeight);
				}
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
				if (weighted) {
					edgeweight w = this->eweights[u][vi];
					handle(u, v, w);
				} else {
					handle(u, v, defaultEdgeWeight);
				}
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
				if (weighted) {
					edgeweight w = this->eweights[u][vi];
					handle(u, v, w);
				} else {
					handle(u, v, defaultEdgeWeight);
				}
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forEdgesWithAttribute_double(int attrId, L handle) {
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
inline void NetworKit::Graph::forEdgesWithAttribute_double(int attrId, L handle) const {
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






/** OLD STUFF **/


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
	if (weighted) {
		for (index i = 0; i < (index) adja[u].size(); ++i) {
			node v = adja[u][i];
			if (v != none) {
				edgeweight ew = eweights[u][i];
				handle(v, ew);
				assert(ew == weight(u, v));
			}
		}
	} else {
		for (index i = 0; i < (index) adja[u].size(); ++i) {
			node v = adja[u][i];
			if (v != none) {
				handle(v, defaultEdgeWeight);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forWeightedNeighborsOf(node u, L handle) const {
	if (weighted) {
		for (index i = 0; i < (index) adja[u].size(); ++i) {
			node v = adja[u][i];
			if (v != none) {
				edgeweight ew = eweights[u][i];
				handle(v, ew);
				assert(ew == weight(u, v));
			}
		}
	} else {
		for (index i = 0; i < (index) adja[u].size(); ++i) {
			node v = adja[u][i];
			if (v != none) {
				handle(v, defaultEdgeWeight);
			}
		}
	}
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

	_Vector<node> neighbors = adja[u];
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

	_Vector<node> neighbors = adja[u];
	std::sort(neighbors.begin(), neighbors.end(), hasSmallerDegree);

	for (node v : neighbors) {
		if (v != none) {
			handle(u, v);
		}
	}
}


template<typename L>
inline void NetworKit::Graph::breadthFirstEdgesFrom(node r, L handle) {
	// TODO: implement BFS iterator for edges
	throw std::runtime_error("TODO");
}


template<typename L>
inline void NetworKit::Graph::forWeightedEdgesOf(node u, L handle) {
	const count asize = (count) adja[u].size();
	for (index i = 0; i < asize; ++i) {
		node v = adja[u][i];
		if (v != none) {
			if (weighted) {
				edgeweight w = this->eweights[u][i];
				handle(u, v, w);
			} else {
				handle(u, v, defaultEdgeWeight);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forWeightedEdgesOf(node u, L handle) const {
	for (index i = 0; i < adja[u].size(); ++i) {
		node v = adja[u][i];
		if (v != none) {
			if (weighted) {
				edgeweight w = this->eweights[u][i];
				handle(u, v, w);
			} else {
				handle(u, v, defaultEdgeWeight);
			}
		}
	}
}

template<typename L>
void NetworKit::Graph::BFSfrom(node r, L handle) {
	std::vector<bool> marked(z);
	std::queue<node> q;
	q.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = q.front();
		q.pop();
		// apply function
		handle(u);
		this->forNeighborsOf(u, [&](node v) {
			if (!marked[v]) {
				q.push(v);
				marked[v] = true;
			}
		});
	} while (!q.empty());
};

template<typename L>
void NetworKit::Graph::DFSfrom(node r, L handle) {
	std::vector<bool> marked(z);
	std::stack<node> stack;
	stack.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = stack.top();
		stack.pop();
		// apply function
		handle(u);
		this->forNeighborsOf(u, [&](node v) {
			if (!marked[v]) {
				stack.push(v);
				marked[v] = true;
			}
		});
	} while (!stack.empty());
};

#endif /* GRAPH_H_ */
