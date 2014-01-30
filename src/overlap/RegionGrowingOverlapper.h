/*
 * RegionGrowingOverlapper.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef REGIONGROWINGOVERLAPPER_H_
#define REGIONGROWINGOVERLAPPER_H_

#include <vector>

#include "Overlapper.h"
#include "../clustering/Clustering.h"

namespace NetworKit {

/**
 * Determines the overlap of multiple partitions by region growing (BFS).
 */
class RegionGrowingOverlapper: public NetworKit::Overlapper {

public:

	RegionGrowingOverlapper();

	virtual ~RegionGrowingOverlapper();

	virtual Clustering run(Graph& G, std::vector<Clustering>& clusterings);

};

} /* namespace NetworKit */
#endif /* REGIONGROWINGOVERLAPPER_H_ */
