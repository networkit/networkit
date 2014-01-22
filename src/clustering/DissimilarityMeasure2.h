/*
 * DissimilarityMeasure2.h
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef DISSIMILARITYMEASURE2_H_
#define DISSIMILARITYMEASURE2_H_

#include "../structures/Partition.h"

namespace NetworKit {


/**
 * Base class for all clustering dissimilarity measures.
 */
class DissimilarityMeasure2 {

public:

	DissimilarityMeasure2();

	virtual ~DissimilarityMeasure2();


	virtual double getDissimilarity(Graph& G, Partition& first, Partition& second) = 0;
};

} /* namespace NetworKit */
#endif /* DISSIMILARITYMEASURE2_H_ */
