/*
 * Modularity.h
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef MODULARITY_H_
#define MODULARITY_H_

#include "QualityMeasure.h"

namespace NetworKit {


/**
 * @ingroup community
 * Modularity is a quality index for community detection. It assigns a quality value in [-0.5, 1.0] to  
 * a partition of a graph which is higher for more modular networks and partitions which better capture 
 * the modular structure.
 * 
 * Modularity is defined as:
 *
 * 	$$mod(\zeta) := \frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)}
 * 	- \frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }$$
 */
class Modularity: public NetworKit::QualityMeasure {
protected:
	double gTotalEdgeWeight;

public:

	/** Default constructor */
	Modularity();

	/**
	 * Returns the Modularity of the given clustering with respect to the graph @a G.
	 *
	 * @param zeta The clustering.
	 * @param G The graph.
	 * @return The modularity.
	 */
	virtual double getQuality(const Partition& zeta, const Graph& G);

	/**
	 * @param totalEdgeWeight Sum of all edge weights in @a G. If specified, it does not
	 *        have to be computed.
	 */
	void setTotalEdgeWeight(double totalEdgeWeight);

};

} /* namespace NetworKit */
#endif /* MODULARITY_H_ */
