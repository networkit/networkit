/*
 * KPathCentrality.h
 *
 *  Created on: 05.10.2014
 *      Author: nemes
 */

#ifndef KPATHCENTRALITY_H_
#define KPATHCENTRALITY_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class KPathCentrality: public NetworKit::Centrality {
public:

	/*
	 * the maximum length of paths
	 * default value ln(n+m)
	 */
	count k;
	/*
	 * value in interval [-0.5, 0.5]
	 * tradeoff between runtime and precision
	 * -0.5: maximum precision, maximum runtime
	 *  0.5: lowest precision, lowest runtime
	 * default value 0.2
	 */
	double alpha;

	/**
	 * Constructs the K-Path Centrality class for the given Graph @a G.
	 *
	 * @param G The graph.
	 * @param alpha tradeoff between precision and runtime.
	 * @param k maximum length of paths.
	 */
	KPathCentrality(const Graph& G, double alpha=0.2, count k=0);

	/**
	* Compute K-Path scores sequential.
	*
	*/
	void run() override;

};

} /* namespace NetworKit */

#endif /* KPATHCENTRALITY_H_ */
