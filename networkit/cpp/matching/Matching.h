/*
 * Matching.h
 *
 *  Created on: 03.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef MATCHING_H_
#define MATCHING_H_

# include "../auxiliary/Log.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup matching
 * FIXME: Could be better to store a reference to the according graph;
 */
class Matching {


public:

	/**
	 * Construct new Matching.
	 *
	 * @param[in]	n 	Maximum number of nodes.
	 */
	Matching(uint64_t n);

	/** Default destructor */
	virtual ~Matching() = default;


	/**
	 * Set two nodes @a u and @a v as each others matching partners.
	 *
	 * @param[in] u node.
	 * @param[in] v node.
	 */
	void match(const node& u, const node& v);


	/**
	 * Reset the two nodes @a u and @a v to unmatched.
	 *
	 * @param[in] u node.
	 * @param[in] v node.
	 */
	void unmatch(const node& u, const node& v);


	/**
	 * Check if node is matched.
	 *
	 * @param[in]	u 	node.
	 * @return @c true if u is matched.
	 */
	bool isMatched(const node& u) const;


	/**
	 * Check if the two nodes @a u and @a v are matched.
	 *
	 * @param[in] u node.
	 * @param[in] v node.
	 */
	bool areMatched(const node& u, const node& v) const;

	/**
	 * Check whether this is a proper matching
	 * in the graph, i.e. no two matched edges are adjacent.
	 *
	 * @paramt[in]	G	A graph.
	 * @param[out]		@c true if this is a proper matching.
	 */
	bool isProper(Graph& G) const;


	/**
	 * Get the number of edges in this matching.
	 * @return Number of edges in matching.
	 */
	count size() const;

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
	edgeweight weight(const Graph& g) const;

protected:

	std::vector<node> data; //!< storage of matching nodes
	count n; //!< number of nodes
};

} /* namespace NetworKit */
#endif /* MATCHING_H_ */
