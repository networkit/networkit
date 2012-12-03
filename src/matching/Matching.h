/*
 * Matching.h
 *
 *  Created on: 03.12.2012
 *      Author: cls
 */

#ifndef MATCHING_H_
#define MATCHING_H_


#include "../graph/Graph.h"

namespace EnsembleClustering {

class Matching {

protected:

	node* array;	//!< array of matching partners

public:

	/**
	 * Construct new matching.
	 *
	 * @param[in]	n 	maximum number of nodes
	 */
	Matching(int n);

	/**
	 * Destructor.
	 */
	virtual ~Matching();


	/**
	 * Set two nodes as eachothers matching
	 * partners.
	 *
	 * @param[in]	u	a node
	 * @param[in]	v	a node
	 *
	 */
	void match(const node& u, const node& v);


	/**
	 *  Index operator.
	 *
	 *  @param[in]	u	a node
	 */
	node& operator[](const node& u);

	/**
	 * Index operator for const instances of this class.
	 *
	 * @param[in]	u 	a node
	 */
	const node& operator[](const node& u) const;


	/**
	 * Check if node is matched
	 *
	 * @param[in]	u 	a node
	 * @param[out]		true if u is matched
	 */
	bool isMatched(const node& u) const;

	/**
	 * Check whether this is a proper matching
	 * in the graph, i.e. no two edges are adjacent
	 *
	 * @paramt[in]	G	a graph
	 * @param[out]		true if this is a proper matching
	 */
	bool isProper(Graph& G) const;


};

} /* namespace EnsembleClustering */
#endif /* MATCHING_H_ */
