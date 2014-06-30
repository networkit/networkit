/*
 * ModularitySequential.h
 *
 *  Created on: 14.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef MODULARITYSEQUENTIAL_H_
#define MODULARITYSEQUENTIAL_H_

#include "QualityMeasure.h"

namespace NetworKit {

/**
 * The ModularitySequential class is the sequential version of Modularity.
 */
class ModularitySequential: public NetworKit::QualityMeasure {
public:

	/** Default destructor */
	virtual ~ModularitySequential();

	/**
	 * Returns the Modularity of the given clustering with respect to the graph @a G.
	 *
	 * @param zeta The clustering.
	 * @param G The graph.
	 * @return The modularity.
	 */
	virtual double getQuality(const Partition& zeta, const Graph& G);

};

} /* namespace NetworKit */
#endif /* MODULARITYSEQUENTIAL_H_ */
