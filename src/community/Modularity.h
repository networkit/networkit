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



public:

	Modularity();

	virtual ~Modularity();

	/**
	 * Returns the Modularity of the given clustering with respect to the graph G.
	 *
	 */
	virtual double getQuality(const Partition& zeta, const Graph& G);
};

} /* namespace NetworKit */
#endif /* MODULARITY_H_ */
