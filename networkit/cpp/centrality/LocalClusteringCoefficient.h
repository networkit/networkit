/*
 * LocalClusteringCoefficient.h
 *
 *  Created on: 31.03.2015
 *      Author: maxv
 */

#ifndef LOCALCLUSTERINGCOEFFICIENT_H_
#define LOCALCLUSTERINGCOEFFICIENT_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class LocalClusteringCoefficient: public NetworKit::Centrality {
public:
	/**
	 * Constructs the LocalClusteringCoefficient class for the given Graph @a G. If the local clustering coefficient scores should be normalized,
	 * then set @a normalized to <code>true</code>.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>true</code> if scores should be normalized in the interval [0,1].
	 * @param computeEdgeCentrality Set this parameter to <code>true</code> if edge betweenness should be computed as well.
	 */
	LocalClusteringCoefficient(const Graph& G, bool normalized=false, bool computeEdgeCentrality=false);



	/**
	* Compute the local clustering coefficient.
	*
	*/
	void run() override;

};

} /* namespace NetworKit */

#endif /* LOCALCLUSTERINGCOEFFICIENT_H_ */
