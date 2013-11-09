/*
 * SampledRandMeasure.h
 *
 *  Created on: 01.07.2013
 *      Author: cls
 */

#ifndef SAMPLEDRANDMEASURE_H_
#define SAMPLEDRANDMEASURE_H_

#include "DissimilarityMeasure.h"

namespace NetworKit {

/**
 * The node-structural Rand measure assigns a similarity value in [0,1]
 * to two partitions of a graph, by considering pairs of nodes.
 * This implementation approximates the index by sampling.
 */
class SampledRandMeasure: public NetworKit::DissimilarityMeasure {

public:

	SampledRandMeasure(count maxSamples);

	virtual ~SampledRandMeasure();

	virtual double getDissimilarity(Graph& G, Clustering& first, Clustering& second);

protected:

	count maxSamples;

};

} /* namespace NetworKit */
#endif /* SAMPLEDRANDMEASURE_H_ */
