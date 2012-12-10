/*
 * Modularity.h
 *
 *  Created on: 10.12.2012
 *      Author: cls
 */

#ifndef MODULARITY_H_
#define MODULARITY_H_

#include <unordered_map>

#include "QualityMeasure.h"

#include "../aux/IndexMap.h"
#include "../graph/Graph.h"
#include "../graph/NodeMap.h"


namespace EnsembleClustering {

class Modularity: public EnsembleClustering::QualityMeasure {

protected:

	NodeMap<double>* incidentWeight;	//<! node -> sum of the weight of incident edges

	/**
	 * Precompute some values depending on the graph instance
	 * to be used in getQuality.
	 */
	virtual void precompute();

public:

	Modularity(Graph& G);

	virtual ~Modularity();

	/**
	 * Returns the Modularity of the given clustering with respect to the graph instance.
	 *
	 * Modularity is defined as:
	 *
	 * 	$$mod(\zeta) := \frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)}
	 * 	- \frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }$$
	 */
	virtual double getQuality(Clustering& zeta);
};

} /* namespace EnsembleClustering */
#endif /* MODULARITY_H_ */
