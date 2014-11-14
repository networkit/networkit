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
class DegreeCentrality: public NetworKit::Centrality {
public:
	/**
	 * Constructs the DegreeCentrality class for the given Graph @a G. If the betweenness scores should be normalized,
	 * then set @a normalized to <code>true</code>.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>true</code> if scores should be normalized in the interval [0,1].
	 */
	DegreeCentrality(const Graph& G, bool normalized=false);

	void run() override;

	/**
	 * @return the theoretical maximum degree centrality, which is $n$ (including the possibility of a self-loop)
	 */
	double maximum();
};

} /* namespace NetworKit */

#endif /* DEGREECENTRALITY_H_ */
