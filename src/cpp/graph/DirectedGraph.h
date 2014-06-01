/*
 * AbstractGraph.cpp
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef DIRECETD_GRAPH_H_
#define DIRECETD_GRAPH_H_

#include <vector>
#include <algorithm>

#include "IDGraph.h"
#include "AbstractGraph.h"
#include "Coordinates.h"
#include "../viz/Point.h"

namespace NetworKit {

/**
 * An directed graph (with optional weights) and parallel iterator methods.
 */
// TODO: add final to class
class DirectedGraph : public IDGraph, public AbstractGraph {

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

	// per node data
	std::vector< NodeDegree > deg; //!< degree of each node (size of neighborhood)

	// per edge data
	std::vector<std::pair<std::vector<node>, std::vector<node> > > adja;
	//!< neighbors/adjacencies, adja.first is adjacencyarray with outgoing edges, adja.second holds incoming edges

	std::vector<std::vector<edgeweight> > eweights; //!< edge weights

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

	/* move assignments and more currently removed because a problems with gcc 4.7.1 on phipute1 */
	// DirectedGraph(const DirectedGraph& other) = default;
	// DirectedGraph(DirectedGraph&& other) = default;
	// DirectedGraph& operator=(DirectedGraph&& other) = default;
	// DirectedGraph& operator=(const DirectedGraph& other) = default;

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
	 * Return the number of neighbors for node v.
	 */
	count degreeIn(node v) const;

	/**
	 * Return the number of outgoing edges from node v.
	 */
	count degreeOut(node v) const;

	/**
	 * @return Weighted degree of @a v.
	 */
	edgeweight weightedDegreeIn(node v) const;

	/**
	 * @return Weighted degree of @a v.
	 */
	edgeweight weightedDegreeOut(node v) const;


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
	// template<typename L> void forEdges(L handle) const;
	void forEdges(FEdge f) const;

	/**
	 * Iterate in parallel over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForEdges(L handle) const;

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
	 */
	template<typename L> void parallelForWeightedEdges(L handle) const;

	/**
	 * Iterate over all edges of the const graph and call handler (lambda closure).
	 *
	 *	@param[in]	attrId		attribute id
	 *	@param[in]	handle 		takes arguments (u, v, a) where a is an edge attribute of edge {u, v}
	 *
	 */
	template<typename L> void forEdgesWithAttribute_double(int attrId, L handle) const;

	
	/** NEIGHBORHOOD ITERATORS **/

	/**
	 * Iterate over all adjacent nodes, which have an edge from u.
	 */
	// template<typename L> void forOutNeighborsOf(node u, L handle) const;
	void forOutNeighborsOf(node u, FNode handle) const;

	/**
	 * Iterate over all adjacent nodes, which have an edge to u.
	 */
	// template<typename L> void forInNeighborsOf(node u, L handle) const;
	void forInNeighborsOf(node u, FNode handle) const;

	/**
	 * Iterate over all outgoing edge weights of a node and call handler (lamdba closure).
	 */
	template<typename L> void forWeightedOutNeighborsOf(node u, L handle) const;
	
	/**
	 * Iterate over all incoming edge weights of a node and call handler (lamdba closure).
	 */
	template<typename L> void forWeightedInNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all outgoing edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forOutEdgesOf(node u, L handle) const;
	
	/**
	 * Iterate over all incoming edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forInEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all incident edges from u and call handler (lamdba closure).
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 *
	 */
	template<typename L> void forWeightedOutEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all incident edges from u and call handler (lamdba closure).
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 *
	 */
	template<typename L> void forWeightedInEdgesOf(node u, L handle) const;


	/** REDUCTION ITERATORS **/

	/**
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForWeightedEdges(L handle) const;

};

} /* namespace NetworKit */


/** EDGE ITERATORS **/

// template<typename L>
// inline void NetworKit::DirectedGraph::forEdges(L handle) const {
inline void NetworKit::DirectedGraph::forEdges(FEdge f) const {
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < (adja[u].first).size(); i++) {
			node v = this->adja[u].first[i];
			if (v != none) {
				f(u, v);
			}
		}
	}
}


template<typename L>
inline void NetworKit::DirectedGraph::parallelForEdges(L handle) const {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < (adja[u].first).size(); i++) {
			node v = this->adja[u].first[i];
			if (v != none) {
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forWeightedEdges(L handle) const {
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < (adja[u].first).size(); i++) {
			node v = this->adja[u].first[i];
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
		for (index i = 0; i < (adja[u].first).size(); i++) {
			node v = this->adja[u].first[i];
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
inline void NetworKit::DirectedGraph::forEdgesWithAttribute_double(int attrId, L handle) const {
	std::vector<std::vector<double> > edgeMap = this->edgeMaps_double[attrId];
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < (adja[u].first).size(); i++) {
			node v = this->adja[u].first[i];
			double attr = edgeMap[u][i];
			if (v != none) {
				handle(u, v, attr);
			}
		}
	}
}


/** NEIGHBORHOOD ITERATORS **/

//template<typename L>
// inline void NetworKit::DirectedGraph::forOutNeighborsOf(node u, L handle) const {
inline void NetworKit::DirectedGraph::forOutNeighborsOf(node u, NetworKit::FNode handle) const {
	for (index i = 0; i < (adja[u].first).size(); i++) {
		node v = this->adja[u].first[i];
		if (v != none) {
			handle(v);
		}
	}
}

 //template<typename L>
// inline void NetworKit::DirectedGraph::forInNeighborsOf(node u, L handle) const {
inline void NetworKit::DirectedGraph::forInNeighborsOf(node u, NetworKit::FNode handle) const {
	for (index i = 0; i < (adja[u].second).size(); i++) {
		node v = this->adja[u].second[i];
		if (v != none) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forWeightedOutNeighborsOf(node u, L handle) const {
	if (weighted) {
		for (index i = 0; i < (adja[u].first).size(); i++) {
			node v = adja[u].first[i];
			if (v != none) {
				edgeweight ew = eweights[u][i];
				handle(v, ew);
				assert(ew == weight(u, v));
			}
		}
	} else {
		for (index i = 0; i < (adja[u].first).size(); i++) {
			node v = adja[u].first[i];
			if (v != none) {
				handle(v, defaultEdgeWeight);
			}
		}
	}
}
	
template<typename L>
inline void NetworKit::DirectedGraph::forWeightedInNeighborsOf(node u, L handle) const {
	if (weighted) {
		for (index i = 0; i < (adja[u].second).size(); i++) {
			node v = adja[u].second[i];
			if (v != none) {
				edgeweight ew = eweights[u][i];
				handle(v, ew);
				assert(ew == weight(u, v));
			}
		}
	} else {
		for (index i = 0; i < (adja[u].second).size(); i++) {
			node v = adja[u].second[i];
			if (v != none) {
				handle(v, defaultEdgeWeight);
			}
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forOutEdgesOf(node u, L handle) const {
	for(index i = 0; i < (adja[u].first).size(); i++) {
		node v = this->adja[u].first[i];
		if (v != none) {
			handle(u, v);
		}
	}
}
	
template<typename L>
inline void NetworKit::DirectedGraph::forInEdgesOf(node u, L handle) const {
	for(index i = 0; i < (adja[u].second).size(); i++) {
		node v = this->adja[u].second[i];
		if (v != none) {
			handle(v, u);
		}
	}
}

template<typename L>
inline void NetworKit::DirectedGraph::forWeightedOutEdgesOf(node u, L handle) const {
	for (index i = 0; i < (adja[u].first).size(); i++) {
		node v = adja[u].first[i];
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
inline void NetworKit::DirectedGraph::forWeightedInEdgesOf(node u, L handle) const {
		for (index i = 0; i < adja[u].second.size(); i++) {
		node v = adja[u].second[i];
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


/** REDUCTION ITERATORS **/

template<typename L>
double NetworKit::DirectedGraph::parallelSumForWeightedEdges(L handle) const {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node u = 0; u < z; u++) {
		for (index i = 0; i < (adja[u].first).size(); i++) {
			node v = this->adja[u].first[i];
			edgeweight ew = this->eweights[u][i];
			if (v != none) {
				sum += handle(u, v, ew);
			}
		}
	}
	return sum;
}

#endif /* DIRECETD_GRAPH_H_ */
