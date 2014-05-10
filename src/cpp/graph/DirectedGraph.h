/*
 * AbstractGraph.cpp
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef DIRECETD_GRAPH_H_
#define DIRECETD_GRAPH_H_

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
 * An directed graph (with optional weights) and parallel iterator methods.
 */
class DirectedGraph final : public AbstractGraph {

protected:

	struct NodeDegree {
		count in;
		count out;
		NodeDegree():
			in(0),
			out(0)
		{}
		count total() const { return in + out; }
	};

	// TODO: structure for Edge (start, end, weight)?

	// per node data
	_Vector< NodeDegree > deg; //!< degree of each node (size of neighborhood)

	// per edge data
	_Vector<_Vector<node> > adja; //!< neighbors/adjacencies, starting with all incoming edges, inOut marks the first outgoing edge
	_Vector<index> inOut; //!< index of first outgoing edge in adja
	_Vector<_Vector<edgeweight> > eweights; //!< edge weights

	// user-defined edge attributes

	/**
	 * Return the index of v in the incoming edges adjacency array of u.
	 */
	index findIn(node u, node v) const;

	/**
	 * Return the index of v in the outgoing edges adjacency array of u.
	 */
	index findOut(node u, node v) const;

public:

	/** 
	 * Create a graph of n nodes.
	 */
	DirectedGraph(count n = 0, bool weighted = false);

	DirectedGraph(const DirectedGraph& other) = default;

	DirectedGraph(DirectedGraph&& other) = default;

	virtual ~DirectedGraph();

	DirectedGraph& operator=(DirectedGraph&& other) = default;

	DirectedGraph& operator=(const DirectedGraph& other) = default;


	/** Only to be used from Cython */
	void stealFrom(DirectedGraph& input);


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
	 * Return the number of neighbors for node v. For directed graphs this is the sum of
	 * in- and outgoing edges.
	 */
	virtual count degree(node v) const;

	/**
	 * Return the number of incoming edges to node v. For an undirected graph this is the
	 * same as degree(v).
	 */
	count degreeIn(node v) const;

	/**
	 * Return the number of outgoing edges from node v. For an undirected graph this is the
	 * same as degree(v).
	 */
	count degreeOut(node v) const;


	/** EDGE MODIFIERS **/

	/**
	 * Insert an directed edge between from @a u to @a v.
	 */
	void addEdge(node u, node v, edgeweight weight = defaultEdgeWeight);

	/**
	 * Check if directed edge {u,v} exists.
	 *
	 */
	bool hasEdge(node u, node v) const;

	/**
	 * Remove directed edge between from @a u to @a v.
	 */
	// void removeEdge(node u, node v);

	/**
	 * Merges edge {u,v} to become a supernode. Edges to u and v are
	 * rewired, multiple edges merged and their weights added.
	 * The vertex weights of @a u and @a v are added.
	 * A self-loop is only created if @a discardSelfLoop is set to false.
	 *
	 * @return New node that has been created if u != v. Otherwise none.
	 */
	// node mergeEdge(node u, node v, bool discardSelfLoop = true);

	/** GLOBAL PROPERTIES **/

	/** 
	 * Return true if this graph supports directed edges.
	 */
	virtual bool isDirected() const;


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

};

} /* namespace NetworKit */


/** EDGE ITERATORS **/

template<typename L>
inline void NetworKit::DirectedGraph::forEdges(L handle) {
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->inOut[u]; i++) {
			node v = this->adja[u][i];
			if (v != none) {
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forEdges(L handle) const {
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->inOut[u]; i++) {
			node v = this->adja[u][i];
			if (v != none) {
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::parallelForEdges(L handle) {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->inOut[u]; i++) {
			node v = this->adja[u][i];
			if (v != none) {
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::parallelForEdges(L handle) const {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->inOut[u]; i++) {
			node v = this->adja[u][i];
			if (v != none) {
				handle(u, v);
			}
		}
	}
}


template<typename L>
inline void NetworKit::DirectedGraph::forWeightedEdges(L handle) {
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->inOut[u]; i++) {
			node v = this->adja[u][i];
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
}

template<typename L>
inline void NetworKit::DirectedGraph::forWeightedEdges(L handle) const {
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->inOut[u]; i++) {
			node v = this->adja[u][i];
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
}

template<typename L>
inline void NetworKit::DirectedGraph::parallelForWeightedEdges(L handle) {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->inOut[u]; i++) {
			node v = this->adja[u][i];
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
}

template<typename L>
inline void NetworKit::DirectedGraph::parallelForWeightedEdges(L handle) const {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->inOut[u]; i++) {
			node v = this->adja[u][i];
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
}


#endif /* DIRECETD_GRAPH_H_ */
