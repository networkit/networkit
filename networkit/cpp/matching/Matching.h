/*
 * Matching.h
 *
 *  Created on: 03.12.2012
 */

#ifndef MATCHING_H_
#define MATCHING_H_

# include "../auxiliary/Log.h"
#include "../graph/Graph.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup matching
 */
class Matching {


public:

	/**
	 * Construct new Matching.
	 *
	 * @param[in]	z	Maximum number of nodes.
	 */
	Matching(count z=0);


	/**
	 * Set two nodes @a u and @a v as each others matching partners.
	 *
	 * @param[in] u node.
	 * @param[in] v node.
	 */
	void match(node u, node v);


	/**
	 * Reset the two nodes @a u and @a v to unmatched.
	 *
	 * @param[in] u node.
	 * @param[in] v node.
	 */
	void unmatch(node u, node v);


	/**
	 * Check if node is matched.
	 *
	 * @param[in]	u 	node.
	 * @return @c true if u is matched.
	 */
	bool isMatched(node u) const;


	/**
	 * Check if the two nodes @a u and @a v are matched together.
	 *
	 * @param[in] u node.
	 * @param[in] v node.
	 */
	bool areMatched(node u, node v) const;

	/**
	 * Check whether this is a proper matching
	 * in the graph, i.e. no two matched edges are adjacent.
	 *
	 * @paramt[in]	G	A graph.
	 * @param[out]		@c true if this is a proper matching.
	 */
	bool isProper(const Graph& G) const;


	/**
	 * Get the number of edges in this matching.
	 * @return Number of edges in matching.
	 */
	count size(const Graph& G) const;

	/**
	 * Get the matched neighbor of @a v if it exists, otherwise @c none.
	 *
	 * @param[in] v node.
	 * @return Mate of @a v if it exists, otherwise none.
	 */
	index mate(node v) const;

	/**
	 * Get total weight of edges in this matching.
	 * @param[in] g The corresponding graph.
	 * @return Total weight of edges in this matching.
	 */
	edgeweight weight(const Graph& G) const;

	/**
	 * Convert matching to a Partition object where matched nodes
	 * belong to the same subset and unmatched nodes belong to a singleton subset.
	 * @return Partition
	 */
	Partition toPartition(const Graph& G) const;

	/**
	 * Get the actual vector storing the data.
	 * @return vector
	 */
	std::vector<node> getVector() const;

protected:

//	const Graph& G;		// reference to graph
	std::vector<node> data; //!< storage of matching nodes
	// count n; //!< number of nodes
};

} /* namespace NetworKit */
#endif /* MATCHING_H_ */
