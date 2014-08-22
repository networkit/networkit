/*
 * Betweenness2.h
 *
 *  Created on: 19.02.2014
 *      Author: cls, ebergamini
 */

#ifndef BETWEENNESS2_H_
#define BETWEENNESS2_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class Betweenness2: public NetworKit::Centrality {
public:
	/**
	 * Constructs the Betweenness class for the given Graph @a G. If the betweenness scores should be normalized,
	 * then set @a normalized to <code>true</code>.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>true</code> if scores should be normalized in the interval [0,1].
	 */
	Betweenness2(const Graph& G, bool normalized=false);



	/**
	* Compute betweenness scores sequential or parallel depending on @a runUnweightedInParallel.
	*
	*/
	void run() override;

};

} /* namespace NetworKit */

#endif /* BETWEENNESS2_H_ */
