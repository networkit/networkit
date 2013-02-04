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




namespace EnsembleClustering {

/** Typedefs **/

typedef int64_t index;
typedef int64_t node;
typedef float edgeweight;

/** Constants **/

int64_t none = -1;


class Graph2 {



protected:

	// per edge data
	std::vector<std::vector<node> > adja; 	//!< neighbors/adjacencies
	std::vector<std::vector<edgeweight> > eweights;	//!< edge weights

	// per node data
	std::vector<int64_t> next;	//!< index of next free slot in adjacency array
	std::vector<int64_t> deg;	//!< degree of each node

	node maxn;	//!< maximum node id / upper bound of node range
	int64_t n;	//!< number of nodes

	// defaults
	edgeweight defaultWeight;

	int64_t find(node u, node v) const;

public:

	Graph2();

	virtual ~Graph2();

	void insertEdge(node u, node v);


	/**
	 * Check if undirected edge {u,v} exists in G
	 *
	 */
	bool hasEdge(node u, node v) const;

	/**
	 * Return edge weight.
	 *
	 * Return 0 if edge does not exist.
	 */
	edgeweight weight(node u, node v) const; // TODO: necessary to inline?

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
	 * Return the number of (non-isolated) nodes in the graph.
	 *
	 * TODO: Maybe this should be changed to support isolated nodes.
	 */
	int64_t numberOfNodes() const;

	/**
	 * Return the number of edges in the graph.
	 */
	int64_t numberOfEdges() const;





	// TODO: lambda iterators

	/**
	 * Iterate over all nodes of the graph and execute callback handletion (lambda closure).
	 */
	template<typename L> void forallNodes(L handle, std::string par="");





	/**
	 * Iterate over all neighbors of a node and execute callback handletion (lamdba closure).
	 */
	template<typename L> void forallNeighborsOf(node v, L handle);




	// TODO: const iterators
};

} /* namespace EnsembleClustering */

template<typename L>
inline void EnsembleClustering::Graph2::forallNeighborsOf(node v, L handle) {

}

template<typename L>
inline void EnsembleClustering::Graph2::forallNodes(L handle,
		std::string par) {
	// TODO: this assumes that no nodes are deleted
	assert ((par == "") || (par == "parallel"));
	#pragma omp parallel for if (par == "parallel")
	for (node v = 1; v <= n; ++v) {
		// call here
		handle(v);
	}
}

#endif /* GRAPH2_H_ */
