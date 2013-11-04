/*
 * ConcurrentGraph.h
 *
 *  Created on: 24.10.2013
 *      Author: cls
 */

#ifndef CONCURRENTGRAPH_H_
#define CONCURRENTGRAPH_H_

#include <cassert>
#include <stdexcept>
#include <tbb/concurrent_vector.h>

#define none std::numeric_limits<index>::max()

namespace NetworKit {

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based
typedef double edgeweight; // edge weight type

/**
* Prototype of a graph data structure which is safe for concurrent modification.
*/
class ConcurrentGraph {
public:

	// defaults
	static constexpr double defaultEdgeWeight = 1.00;
	static constexpr edgeweight nullWeight = 0.0;

	ConcurrentGraph();

	ConcurrentGraph(count n);

	ConcurrentGraph(const ConcurrentGraph& other) = default;

	ConcurrentGraph(ConcurrentGraph&& other) = default;

	virtual ~ConcurrentGraph();

	/**
	 * Assignment operator
	 */
	ConcurrentGraph& operator=(ConcurrentGraph&& other) = default;

	/**
	 * Assignment operator
	 */
	ConcurrentGraph& operator=(const ConcurrentGraph& other) = default;


	/**
	 * Return the number of edges in the graph.
	 *
	 *	 */
	count numberOfEdges() const;


	/**
	 * @return Number of neighbors.
	 */
	count degree(node v) const;

	/**
	 * Insert an undirected edge between two nodes.
	 */
	void addEdge(node u, node v, edgeweight weight = defaultEdgeWeight);

	/**
	 * Remove undirected edge between two nodes.
	 */
	void removeEdge(node u, node v);


private:
	// scalars
	count n; //!< current number of nodes
	count m; //!< current number of edges
	node z; //!< current upper bound of node ids

	// per node data
	tbb::concurrent_vector<count> deg; //!< degree of each node (size of neighborhood)

	// per edge data
	tbb::concurrent_vector<tbb::concurrent_vector<node> > adja; //!< neighbors/adjacencies
	tbb::concurrent_vector<tbb::concurrent_vector<node> > eweights; //!< edge weights


	/**
	 * Return the index of v in the adjacency array of u.
	 */
	index find(node u, node v) const;
};

} /* namespace NetworKit */
#endif /* CONCURRENTGRAPH_H_ */
