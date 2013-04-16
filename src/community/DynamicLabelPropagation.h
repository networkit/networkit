/*
 * DynamicLabelPropagation.h
 *
 *  Created on: 27.03.2013
 *      Author: cls
 */

#ifndef DYNAMICLABELPROPAGATION_H_
#define DYNAMICLABELPROPAGATION_H_

#include <algorithm>

#include "DynamicClusterer.h"
#include "../aux/Timer.h"

namespace NetworKit {

typedef cluster label;

class DynamicLabelPropagation: public NetworKit::DynamicClusterer {

protected:

	count n; //!< current number of nodes TOOD: update
	Clustering labels;
	std::vector<bool> activeNodes;
	std::vector<double> weightedDegree; //!< precompute and update weighted degree for performance reasons
	count updateThreshold;
	count nUpdated; //!< number of nodes updated in last iteration (?)

public:

	/**
	 * @param[in]	G		graph
	 * @param[in]	theta	update threshold
	 */
	DynamicLabelPropagation(Graph& G, count theta);

	virtual ~DynamicLabelPropagation();

	virtual Clustering run();

	virtual std::string toString() const;

	// GRAPH DYNAMICS INTERFACE

	virtual void onNodeAddition(node u);

	virtual void onNodeRemoval(node u);

	virtual void onEdgeAddition(node u, node v);

	virtual void onEdgeRemoval(node u, node v);

	virtual void onWeightUpdate(node u, node v, edgeweight w);

};

} /* namespace NetworKit */
#endif /* DYNAMICLABELPROPAGATION_H_ */
