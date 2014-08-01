/*
 * DissimilarityMeasure.h
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef DISSIMILARITYMEASURE_H_
#define DISSIMILARITYMEASURE_H_

#include "../structures/Partition.h"
#include "../structures/Cover.h"

namespace NetworKit {


/**
 * @ingroup community
 * Base class for all clustering dissimilarity measures.
 */
class DissimilarityMeasure {

public:

	virtual double getDissimilarity(const Graph& G, const Partition& first, const Partition& second) = 0;


	virtual double getDissimilarity(const Graph &G, const Cover &first, const Cover &second);
};

} /* namespace NetworKit */
#endif /* DISSIMILARITYMEASURE_H_ */
