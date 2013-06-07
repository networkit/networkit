/*
 * LocalQualityMeasure.
h
 *
 *  Created on: 28.05.2013
 *      Author: Yassine Marrakchi
 */

#ifndef LOCALQUALITYMEASURE_H_
#define LOCALQUALITYMEASURE_H_

#include <unordered_set>

#include "../clustering/Clustering.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Abstract base class for all local clustering quality measures.
 */
class LocalQualityMeasure {


public:

	LocalQualityMeasure();

	virtual ~LocalQualityMeasure();

	virtual double getQuality(std::unordered_set<node>& C, Graph& G) = 0;
};

} /* namespace NetworKit */


#endif /* LOCALQUALITYMEASURE_H_ */
