/*
 * SelectiveDissimilarityMeasure.h
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef SELECTIVEDISSIMILARITYMEASURE_H_
#define SELECTIVEDISSIMILARITYMEASURE_H_

#include <unordered_set>

#include "../clustering/Clustering.h"

namespace NetworKit {

/**
 * Abstract base class for all measures which evaluate a single community
 * with respect to a given ground-truth clustering of the graph.
 */
class SelectiveDissimilarityMeasure {

public:

	SelectiveDissimilarityMeasure();

	virtual ~SelectiveDissimilarityMeasure();

	virtual double getDissimilarity(const std::unordered_set<node>& community, const Clustering& groundTruth) = 0;
};

} /* namespace NetworKit */
#endif /* SELECTIVEDISSIMILARITYMEASURE_H_ */
