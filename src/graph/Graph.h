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

extern "C" {
#include "stinger.h"
}


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

#define EDGE_SOURCE STINGER_EDGE_SOURCE
#define EDGE_DEST STINGER_EDGE_DEST

#define FORALL_EDGES_OF_NODE_BEGIN(G, V) STINGER_FORALL_EDGES_OF_VTX_BEGIN(G.asSTINGER(), V)
#define FORALL_EDGES_OF_NODE_END() STINGER_FORALL_EDGES_OF_VTX_END()

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
	stinger* asSTINGER();

	/**
	 * Insert a weighted, undirected edge.
	 */
	void insertEdge(node u, node v, double weight=defaultEdgeWeight, int64_t type=defaultEdgeType, int64_t timestamp=defaultTimeStamp);

	/**
	 * Check if undirected edge {u,v} exists in G
	 *
	 */
	bool hasEdge(node u, node v);


	/**
	 * Return node weight.
	 */
	double getWeight(node v);

	/**
	 * Return edge weight.
	 */
	double getWeight(edge uv);

	/**
	 * Return edge weight.
	 * Equivalent to getWeight(edge uv)
	 */
	double getWeight(node u, node v);


	/**
	 * Return the degree (number of incident edges).
	 */
	int64_t getDegree(node u);

	/**
	 * Return the number of edges in the graph.
	 */
	int64_t numberOfEdges();

	/**
	 * Return the number of (non-isolated) nodes in the graph.
	 *
	 * TODO: Maybe this should be changed to support isolated nodes.
	 */
	int64_t numberOfNodes();


	/**
	 * Get the first node index (for iteration over all nodes)
	 */
	node firstNode();


	/**
	 * Get the last node index (for iteration over all nodes).
	 */
	node lastNode();



	template<typename Callback> void forallEdges(bool parallel, Callback callback);



};

} /* namespace EnsembleClustering */

template<typename Callback>
inline void EnsembleClustering::Graph::forallEdges(bool parallel, Callback callback) {
	STINGER_FORALL_EDGES_BEGIN(this->stingerG, this->defaultEdgeType) {
		node u = STINGER_EDGE_SOURCE;
		node v = STINGER_EDGE_DEST;
		// call the supplied callback
		callback(u, v);
	} STINGER_FORALL_EDGES_END();

}

#endif /* GRAPH_H_ */
