/*
 * Closeness.h
 *
 *  Created on: 03.10.2014
 *      Author: nemes
 */

#ifndef CLOSENESS_H_
#define CLOSENESS_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class Closeness: public NetworKit::Centrality {
public:
	/**
	 * Constructs the Betweenness class for the given Graph @a G. If the betweenness scores should be normalized,
	 * then set @a normalized to <code>true</code>.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>true</code> if scores should be normalized in the interval [0,1].
	 */
	Closeness(const Graph& G, bool normalized=false);



	/**
	* Compute betweenness scores sequential or parallel depending on @a runUnweightedInParallel.
	*
	*/
	void run() override;

	double centralization();
};

} /* namespace NetworKit */

#endif /* CLOSENESS_H_ */
