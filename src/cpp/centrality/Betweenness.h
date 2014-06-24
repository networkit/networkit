/*
 * Betweenness.h
 *
 *  Created on: 19.02.2014
 *      Author: hm
 */

#ifndef BETWEENNESS_H_
#define BETWEENNESS_H_

#include "Centrality.h"

namespace NetworKit {

class Betweenness: public NetworKit::Centrality {
public:
	/**
	 * Constructs the Betweenness class for the given Graph @a G. If the betweenness scores should be normalized,
	 * then set @a normalized to <code>true</code>.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>true</code> if scores should be normalized in the interval [0,1].
	 */
	Betweenness(const Graph& G, bool normalized=false);

	/**
	 * Compute betweenness scores sequential or parallel depending on @a runUnweightedInParallel.
	 *
	 * @param runUnweightedInParallel If set to <code>true</code> the computation is done in parallel.
	 */
	void run(bool runUnweightedInParallel);

	void run() override { run(false); }

};

} /* namespace NetworKit */

#endif /* BETWEENNESS_H_ */
