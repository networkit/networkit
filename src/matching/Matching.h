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
#include "../graph/NodeMap.h"

namespace NetworKit {

/**
 * FIXME: Could be better to store a reference to the according graph;
 */
class Matching : public NodeMap<node> {


public:

	/**
	 * Construct new matching.
	 *
	 * @param[in]	n 	maximum number of nodes
	 */
	Matching(uint64_t n);

	/**
	 * Destructor.
	 */
	virtual ~Matching();


	/**
	 * Set two nodes as eachothers matching
	 * partners.
	 *
	 *
	 */
	void match(const node& u, const node& v);


	/**
	 * Reset the two nodes to unmatched.
	 */
	void unmatch(const node& u, const node& v);


	/**
	 * Check if node is matched
	 *
	 * @param[in]	u 	a node
	 * @param[out]		true if u is matched
	 */
	bool isMatched(const node& u) const;


	/**
	 * Check if the two nodes are matched.
	 *
	 */
	bool areMatched(const node& u, const node& v) const;

	/**
	 * Check whether this is a proper matching
	 * in the graph, i.e. no two edges are adjacent.
	 *
	 *
	 * @paramt[in]	G	a graph
	 * @param[out]		true if this is a proper matching
	 */
	bool isProper(Graph& G) const;


	/**
	 * @return Number of edges in matching.
	 */
	count size() const;

	/**
	 * @return Mate of @a v if it exists, otherwise none.
	 */
	index mate(node v) const;

	/**
	 * @return Total weight of edges in the matching.
	 */
	edgeweight weight(const Graph& g) const;
};

} /* namespace NetworKit */
#endif /* MATCHING_H_ */
