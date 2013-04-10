/*
 * DynamicClusterer.h
 *
 *  Created on: 27.03.2013
 *      Author: cls
 */

#ifndef DYNAMICCLUSTERER_H_
#define DYNAMICCLUSTERER_H_

#include "Clusterer.h"

namespace NetworKit {

class DynamicClusterer {

public:

	DynamicClusterer(Graph& G);

	virtual ~DynamicClusterer();

	/**
	 * Produce clustering.
	 */
	virtual Clustering run() = 0;

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;


	// GRAPH DYNAMICS INTERFACE

	virtual void onNodeAddition(node u) = 0;

	virtual void onNodeRemoval(node u) = 0;

	virtual void onEdgeAddition(node u, node v) = 0;

	virtual void onEdgeRemoval(node u, node v) = 0;

	virtual void onWeightUpdate(node u, node v, edgeweight w) = 0;


};

} /* namespace NetworKit */
#endif /* DYNAMICCLUSTERER_H_ */
