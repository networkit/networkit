/*
 * Graph.h
 *
 *  Created on: 28.11.2012
 *      Author: cls
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <utility>
#include <cinttypes>
#include <string>
#include <queue>

extern "C" {
#include "stinger.h"
}

#include "../aux/log.h"

namespace EnsembleClustering {


/** Typedefs **/


typedef int64_t node; //!< a node is an integer logical index. it is 1-based!
typedef std::pair<node, node> edge; //!< an undirected edge is a pair of nodes (indices)



/** Traversal macros
 *
 * These are modified versions of the macros defined in stinger/include/stinger-traversal.h
 *
 *
 *  **/


#define FORALL_EDGES_BEGIN(G) STINGER_FORALL_EDGES_BEGIN(G.asSTINGER(), G.defaultEdgeType)
#define FORALL_EDGES_END() STINGER_FORALL_EDGES_END()

#define PARALLEL_FORALL_EDGES_BEGIN(G) STINGER_PARALLEL_FORALL_EDGES_BEGIN(G.asSTINGER(), G.defaultEdgeType)
#define PARALLEL_FORALL_EDGES_END() STINGER_PARALLEL_FORALL_EDGES_END()

#define READ_ONLY_FORALL_EDGES_BEGIN(G) STINGER_READ_ONLY_FORALL_EDGES_BEGIN(G.asSTINGER(), G.defaultEdgeType)
#define READ_ONLY_FORALL_EDGES_END() STINGER_READ_ONLY_FORALL_EDGES_END()

#define READ_ONLY_PARALLEL_FORALL_EDGES_BEGIN(G) STINGER_READ_ONLY_PARALLEL_FORALL_EDGES_BEGIN(G.asSTINGER(), G.defaultEdgeType)
#define READ_ONLY_PARALLEL_FORALL_EDGES_END() STINGER_READ_ONLY_PARALLEL_FORALL_EDGES_END()

#define EDGE_SOURCE STINGER_EDGE_SOURCE
#define EDGE_DEST STINGER_EDGE_DEST

#define FORALL_EDGES_OF_NODE_BEGIN(G, V) STINGER_FORALL_EDGES_OF_VTX_BEGIN(G.asSTINGER(), V)
#define FORALL_EDGES_OF_NODE_END() STINGER_FORALL_EDGES_OF_VTX_END()

#define READ_ONLY_FORALL_EDGES_OF_NODE_BEGIN(G, V) STINGER_READ_ONLY_FORALL_EDGES_OF_VTX_BEGIN(G.asSTINGER(), V)
#define READ_ONLY_FORALL_EDGES_OF_NODE_END() STINGER_READ_ONLY_FORALL_EDGES_OF_VTX_END()

// TODO: add the missing macros



/** Graph interface **/


/**
 * Graph encapsulates a STINGER graph object and provides
 * a more concise interface to it.
 *
 * The graph concept modelled is
 * - undirected
 * - weighted
 * - without self-loops (use node weights instead)
 *
 *
 */
class Graph {

protected:

	stinger* stingerG;
	// TODO: is it necessary to store isolated nodes?
//	int64_t maxNode; 	//!< the maximum index for a node in the graph.
//						//!< Needed to keep track of the number of (possibly isolated) nodes,
//						//!< since stinger does not store isolated nodes


public:
	/** default parameters ***/

	static constexpr double defaultEdgeWeight = 1.0;
	static const int64_t defaultEdgeType = 0;
	static const int64_t defaultTimeStamp = 0;

	/** methods **/


	/**
	 * Construct Graph object with new STINGER graph inside.
	 */
	Graph();

	/**
	 * Initialize with STINGER graph.
	 *
	 * @param[in]	stingerG	a STINGER graph struct
	 */
	Graph(stinger* stingerG);

	~Graph();

	/**
	 * Return the internal STINGER data structure.
	 *
	 */
	stinger* asSTINGER() const;

	/**
	 * Insert a weighted, undirected edge.
	 */
	void insertEdge(node u, node v, double weight=defaultEdgeWeight, int64_t type=defaultEdgeType, int64_t timestamp=defaultTimeStamp);


	/**
	 * Remove an unirected edges.
	 */
	void removeEdge(node u, node v);


	/**
	 * Check if undirected edge {u,v} exists in G
	 *
	 */
	bool hasEdge(node u, node v) const;


	/**
	 * Return node weight.
	 */
	double weight(node v) const;

	/**
	 * Return edge weight.
	 */
	double weight(edge uv) const;

	/**
	 * Return edge weight.
	 *
	 * Return 0 if edge does not exist.
	 */
	inline double weight(node u, node v) const {
		return stinger_edgeweight(this->stingerG, u, v, this->defaultEdgeType);
	}


	/**
	 * Set the weight of an edge
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	inline void setWeight(node u, node v, double weight) {
		// set for both directed edges
		stinger_set_edgeweight(this->stingerG, u, v, this->defaultEdgeType, weight);
		stinger_set_edgeweight(this->stingerG, v, u, this->defaultEdgeType, weight);
	}

	/**
	 * Set the weight of a node.
	 */
	inline void setWeight(node u, double weight) {
		stinger_set_vweight(this->stingerG, u, weight);
	}


	/**
	 * Get the sum of the weight of all edges.
	 */
	double totalEdgeWeight();





	/**
	 * Return the degree (number of incident edges).
	 */
	int64_t degree(node u) const;

	/**
	 * Return the number of edges in the graph.
	 */
	int64_t numberOfEdges() const;


	/**
	 * Add a new node to the graph and return it.
	 */
	node addNode();

	/**
	 * Return the number of (non-isolated) nodes in the graph.
	 *
	 * TODO: Maybe this should be changed to support isolated nodes.
	 */
	int64_t numberOfNodes() const;


	/**
	 * Get the first node index (for iteration over all nodes)
	 */
	node firstNode() const;


	/**
	 * Get the last node index (for iteration over all nodes).
	 */
	node lastNode() const;


	/********** ATTRIBUTES **************/

	// TODO: template <typename T> void addNodeMap(NodeMap<T>& map)



	/********** ITERATION / TRAVERSAL ***********/



	/**
	 * Iterate over all undirected edges of a graph and execute callback function.
	 *
	 * If the callback function is a lambda closure which captures the variables
	 * of the curent scope, iteration looks similar to ordinary loops with STINGER traversal macros.
	 */
	template<typename Callback> void forallEdges(Callback func, std::string par="", std::string write="");


	/**
	 * Iterate over all nodes of the graph and execute callback function (lambda closure).
	 */
	template<typename Callback> void forallNodes(Callback func, std::string par="");


	/**
	 * Iterate over all neighbors of a node and execute callback function (lamdba closure).
	 */
	template<typename Callback> void forallNeighborsOf(node v, Callback func);


	/**
	 * Iterate over all incident edges of a node and execute callback function (lamdba closure).
	 */
	template<typename Callback> void forallEdgesOf(node v, Callback func);


	/**
	 * Iterate over nodes in breadth-first search order starting from r until connected component
	 * of r has been visited.
	 */
	template<typename Callback> void breadthFirstNodesFrom(node r, Callback func);


	/**
	 * Iterate over edges in breadth-first search order starting from node r until connected component
	 * of r has been visited.
	 */
	template<typename Callback> void breadthFirstEdgesFrom(node r, Callback func);


};

} /* namespace EnsembleClustering */

template<typename Callback>
inline void EnsembleClustering::Graph::forallEdges(Callback func, std::string par, std::string write) {
	if (par == "parallel") {
		if (write == "readonly") {
			// parallel, readonly
			STINGER_READ_ONLY_PARALLEL_FORALL_EDGES_BEGIN(this->stingerG, this->defaultEdgeType) {
				node u = STINGER_EDGE_SOURCE;
				node v = STINGER_EDGE_DEST;
				if (u < v) {
					// consider only undirected edges
					// call the supplied callback
					func(u, v);
				}
			} STINGER_READ_ONLY_PARALLEL_FORALL_EDGES_END();
			return;
		} else {
			// parallel, writable
			STINGER_PARALLEL_FORALL_EDGES_BEGIN(this->stingerG, this->defaultEdgeType) {
				node u = STINGER_EDGE_SOURCE;
				node v = STINGER_EDGE_DEST;
				if (u < v) {
					func(u, v);
				}
			} STINGER_PARALLEL_FORALL_EDGES_END();
			return;
		}
	} else {
		if (write == "readonly") {
			// sequential, readonly
			STINGER_READ_ONLY_FORALL_EDGES_BEGIN(this->stingerG, this->defaultEdgeType) {
				node u = STINGER_EDGE_SOURCE;
				node v = STINGER_EDGE_DEST;
				if (u < v) {
					// consider only undirected edges
					func(u, v);
				}
			} STINGER_READ_ONLY_FORALL_EDGES_END();
			return;
		} else {
			// sequential, writable
			STINGER_FORALL_EDGES_BEGIN(this->stingerG, this->defaultEdgeType) {
				node u = STINGER_EDGE_SOURCE;
				node v = STINGER_EDGE_DEST;
				if (u < v) {
					// consider only undirected edges
					func(u, v);
				}
			} STINGER_FORALL_EDGES_END();
			return;
		}
	}


}

template<typename Callback>
inline void EnsembleClustering::Graph::forallNodes(Callback func, std::string par) {
	int64_t n  = this->numberOfNodes();
	#pragma omp parallel for if (par == "parallel")
	for (node v = this->firstNode(); v <= n; ++v) {
		// call node function
		func(v);
	}
}

template<typename Callback>
inline void EnsembleClustering::Graph::forallNeighborsOf(node v, Callback func) {
	STINGER_FORALL_EDGES_OF_VTX_BEGIN(this->stingerG, v) {
		// filtering edges is not necessary because only out-edges are considered by stinger
		node w = EDGE_DEST;
		// call node function
		func(w);
	} STINGER_FORALL_EDGES_OF_VTX_END();
}


template<typename Callback>
inline void EnsembleClustering::Graph::forallEdgesOf(node u, Callback func) {
	STINGER_FORALL_EDGES_OF_VTX_BEGIN(this->stingerG, u) {
		// filtering edges is not necessary because only out-edges are considered by stinger
		node v = STINGER_EDGE_SOURCE; // = u
		node w = STINGER_EDGE_DEST;
		// call edge function
		func(v, w);
	} STINGER_FORALL_EDGES_OF_VTX_END();
}

template<typename Callback>
inline void EnsembleClustering::Graph::breadthFirstNodesFrom(node r, Callback func) {
	std::queue<node> q;
	int64_t n = this->numberOfNodes();
	bool marked[n+1];
	for (int64_t i = 0; i <= n+1; ++i) {
		marked[i] = false;
	}

	q.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = q.front();
		q.pop();
		// apply function
		func(u);
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(this->stingerG, u) {
			// filtering edges is not necessary because only out-edges are considered by stinger
			node v = STINGER_EDGE_DEST;
			if (! marked[v]) {
				DEBUG("pushing node " << v);
				q.push(v);
				marked[v] = true;
			}
		} STINGER_FORALL_EDGES_OF_VTX_END();

	} while (! q.empty());

}

template<typename Callback>
inline void EnsembleClustering::Graph::breadthFirstEdgesFrom(node r, Callback func) {
	std::queue<node> q;
	int64_t n = this->numberOfNodes();
	bool marked[n+1];
	for (int64_t i = 0; i <= n+1; ++i) {
		marked[i] = false;
	}

	DEBUG("pushing node " << r);
	q.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = q.front();
		q.pop();
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(this->stingerG, u) {
			node v = STINGER_EDGE_DEST;
			if (u < v) {
				DEBUG("visiting edge (" << u << "," << v << ")");
				// apply edge function
				func(u, v);
				if (! marked[v]) {
					DEBUG("pushing node " << v);
					q.push(v);
					marked[v] = true;
				}
			}
		} STINGER_FORALL_EDGES_OF_VTX_END();
	} while (! q.empty());
}

#endif /* GRAPH_H_ */
