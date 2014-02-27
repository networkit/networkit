/*
 * SampledRandMeasure.h
 *
 *  Created on: 01.07.2013
 *      Author: cls
 */

#ifndef SAMPLEDGRAPHSTRUCTURALRANDMEASURE_H_
#define SAMPLEDGRAPHTRUCTURALRANDMEASURE_H_

#include "DissimilarityMeasure.h"

namespace NetworKit {

/**
 * The graph-structural Rand measure assigns a similarity value in [0,1]
 * to two partitions of a graph, by considering connected pairs of nodes.
 * This implementation approximates the index by sampling.
 */
class SampledGraphStructuralRandMeasure: public NetworKit::DissimilarityMeasure {

public:

	SampledGraphStructuralRandMeasure(count maxSamples);

	virtual ~SampledGraphStructuralRandMeasure();

	virtual double getDissimilarity(Graph& G, Partition& first, Partition& second);

protected:

	count maxSamples;

};

} /* namespace NetworKit */
#endif /* SAMPLEDRANDMEASURE_H_ */
