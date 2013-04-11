/*
 * DynamicLabelPropagation.h
 *
 *  Created on: 27.03.2013
 *      Author: cls
 */

#ifndef DYNAMICLABELPROPAGATION_H_
#define DYNAMICLABELPROPAGATION_H_

#include "DynamicClusterer.h"

namespace NetworKit {

class DynamicLabelPropagation: public NetworKit::DynamicClusterer {

protected:

	Graph* G;
	std::vector<bool> activeNodes;
	count updateThreshold;
	// Clustering labels;

public:

	DynamicLabelPropagation(Graph& G);

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
