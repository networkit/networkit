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
	 */
	LocalClusteringCoefficient(const Graph& G);



	/**
	* Compute the local clustering coefficient.
	*
	*/
	void run() override;


	/**
	* Get the theoretical maximum of centrality score in the given graph.
	*
	* @return The maximum centrality score.
	*/
	virtual double maximum();

};

} /* namespace NetworKit */

#endif /* LOCALCLUSTERINGCOEFFICIENT_H_ */
