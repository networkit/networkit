/*
 * DegreeCentrality.h
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#ifndef DEGREECENTRALITY_H_
#define DEGREECENTRALITY_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * Node centrality index which ranks nodes by their degree.
 * Optional normalization by maximum degree.
 */
class DegreeCentrality: public Centrality {
public:
	/**
	 * Constructs the DegreeCentrality class for the given Graph @a G. If the centrality scores should be normalized,
	 * then set @a normalized to <code>true</code>. The run() method runs in O(m) time, where m is the number of
	 * edges in the graph.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>true</code> if scores should be normalized in the interval [0,1].
	 * @param outDeg If set to <code>true</code>, computes the centrality based on the out-degrees, otherwise based on the in-degrees.
	 * @param ignoreSelfLoops If set to <code>true</code>, self loops will not be taken into account.
	 */
	DegreeCentrality(const Graph& G, bool normalized=false, bool outDeg=true, bool ignoreSelfLoops=true);

	/**
	 * Computes degree centraity on the graph passed in constructor.
	 */
	void run() override;

	/**
	 * @return the theoretical maximum degree centrality, which is $n$ (including the possibility of a self-loop)
	 */
	double maximum() override;

private:
	bool outDeg, ignoreSelfLoops;
};

} /* namespace NetworKit */

#endif /* DEGREECENTRALITY_H_ */
