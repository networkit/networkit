/*
 * KPath.h
 *
 *  Created on: 05.10.2014
 *      Author: nemes
 */

#ifndef KPATH_H_
#define KPATH_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class KPath: public NetworKit::Centrality {
public:
	/**
	 * Constructs the Betweenness class for the given Graph @a G. If the betweenness scores should be normalized,
	 * then set @a normalized to <code>true</code>.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>true</code> if scores should be normalized in the interval [0,1].
	 */
	KPath(const Graph& G, bool normalized=false);



	/**
	* Compute betweenness scores sequential or parallel depending on @a runUnweightedInParallel.
	*
	*/
	void run() override;

};

} /* namespace NetworKit */

#endif /* KPATH_H_ */
