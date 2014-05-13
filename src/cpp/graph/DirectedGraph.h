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

#include "IDGraph.h"
#include "AbstractGraph.h"
#include "Coordinates.h"
#include "../auxiliary/Log.h"
#include "../Globals.h"
#include "../viz/Point.h"

namespace NetworKit {

/**
 * An directed graph (with optional weights) and parallel iterator methods.
 */
class DirectedGraph final : public IDGraph, public AbstractGraph {

protected:

	// per node data
	_Vector< NodeDegree > deg; //!< degree of each node (size of neighborhood)

	// per edge data
	_Vector<_Vector<node> > adja; //!< neighbors/adjacencies, starting with all incoming edges, inOut marks the first outgoing edge
	_Vector<index> inOut; //!< index of first outgoing edge in adja
	_Vector<_Vector<edgeweight> > eweights; //!< edge weights

	// user-defined edge attributes
	std::vector<std::vector<std::vector<double> > > edgeMaps_double; // contains edge maps (u, v) -> double
	std::vector<double> edgeAttrDefaults_double; // stores default value for edgeMaps_double[i] at index i

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

	virtual count getMemoryUsage() {
		count mem = AbstractGraph::getMemoryUsage();
		mem += sizeof(deg) + sizeof(NodeDegree) * deg.capacity();
		
		mem += sizeof(adja) + sizeof(_Vector<node>) * adja.capacity();
		for (auto& a : adja) {
			mem += sizeof(node) * a.capacity();
		}

		mem += sizeof(inOut) + sizeof(index) * inOut.capacity();

		mem += sizeof(eweights) + sizeof(_Vector<edgeweight>) * eweights.capacity();
		for (auto& w : eweights) {
			mem += sizeof(edgeweight) * w.capacity();
		}

		// TODO attribute stuff

		return mem;
	}

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
	 * @return true if the node is isolated (= degree is 0)
	 */
	bool isIsolated(node v) const;

	/**
	 * Return the number of neighbors for node v.
	 */
	NodeDegree degree(node v) const;

	/** EDGE MODIFIERS **/

	/**
	 * Insert an directed edge between from @a u to @a v.
	 */
	void addEdge(node u, node v, edgeweight weight = defaultEdgeWeight);

	/**
	 * Remove directed edge between from @a u to @a v.
	 */
	void removeEdge(node u, node v);

	/**
	 * Check if directed edge {u,v} exists.
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
	// node mergeEdge(node u, node v, bool discardSelfLoop = true);


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


	/** GLOBAL PROPERTIES **/

	/** 
	 * Return true if this graph supports directed edges.
	 */
	bool isDirected() const;


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

	
	/** ITERATE OVER NEIGHBORS */

	/**
	 * Iterate over all Outgoing edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forOutEdgesOf(node u, L handle);

	/**
	 * Iterate over all Outgoing edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forOutEdgesOf(node u, L handle) const;
	
	/**
	 * Iterate over all Incoming edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forInEdgesOf(node u, L handle);
	
	/**
	 * Iterate over all Incoming edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forInEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all adjacent nodes, which are adjacent to an inedge of u
	 */
	template<typename L> void forOutNeighborsOf(node u, L handle);

	/**
	 * Iterate over all adjacent nodes, which are adjacent to an inedge of u
	 */
	template<typename L> void forOutNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all adjacent nodes, which are adjacent to an outedge of u
	 */
	template<typename L> void forInNeighborsOf(node u, L handle);
	
	/**
	 * Iterate over all adjacent nodes, which are adjacent to an outedge of u
	 */
	template<typename L> void forInNeighborsOf(node u, L handle) const;
	

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

template<typename L>
inline void NetworKit::DirectedGraph::forEdgesWithAttribute_double(int attrId, L handle) {
	std::vector<std::vector<double> > edgeMap = this->edgeMaps_double[attrId];
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->inOut[u]; i++) {
			node v = this->adja[u][i];
			double attr = edgeMap[u][i];
			if (v != none) {
				handle(u, v, attr);
			}
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forEdgesWithAttribute_double(int attrId, L handle) const {
	std::vector<std::vector<double> > edgeMap = this->edgeMaps_double[attrId];
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->inOut[u]; i++) {
			node v = this->adja[u][i];
			double attr = edgeMap[u][i];
			if (v != none) {
				handle(u, v, attr);
			}
		}
	}
}


/** ITERATE OVER NEIGHBORS */

template<typename L>
inline void NetworKit::DirectedGraph::forOutEdgesOf(node u, L handle) {
	for(index i = this->inOut[u]; i < this->adja[u].size(); i++) {
		node v = this->adja[u][i];
		if (v != none) {
			handle(u, v);
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forOutEdgesOf(node u, L handle) const {
	for(index i = this->inOut[u]; i < this->adja[u].size(); i++) {
		node v = this->adja[u][i];
		if (v != none) {
			handle(u, v);
		}
	}
}
	
template<typename L>
inline void NetworKit::DirectedGraph::forInEdgesOf(node u, L handle) {
	for(index i = 0; i < this->inOut[u]; i++) {
		node v = this->adja[u][i];
		if (v != none) {
			handle(u, v);
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forInEdgesOf(node u, L handle) const {
	for(index i = 0; i < this->inOut[u]; i++) {
		node v = this->adja[u][i];
		if (v != none) {
			handle(u, v);
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forOutNeighborsOf(node u, L handle) {
	for (index i = this->inOut[u]; i < this->adja[u].size(); i++) {
		node v = this->adja[u][i];
		if (v != none) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forOutNeighborsOf(node u, L handle) const {
	for (index i = this->inOut[u]; i < this->adja[u].size(); i++) {
		node v = this->adja[u][i];
		if (v != none) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forInNeighborsOf(node u, L handle) {
	for (index i = 0; i < this->inOut[u]; i++) {
		node v = this->adja[u][i];
		if (v != none) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forInNeighborsOf(node u, L handle) const {
	for (index i = 0; i < this->inOut[u]; i++) {
		node v = this->adja[u][i];
		if (v != none) {
			handle(v);
		}
	}
}
	




#endif /* DIRECETD_GRAPH_H_ */
