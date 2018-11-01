/*
 * SampledRandMeasure.h
 *
 *  Created on: 01.07.2013
 *      Author: cls
 */

#ifndef SAMPLEDNODESTRUCTURALRANDMEASURE_H_
#define SAMPLEDNODESTRUCTURALRANDMEASURE_H_

#include "DissimilarityMeasure.h"

namespace NetworKit {

/**
 * @ingroup community
 * The node-structural Rand measure assigns a similarity value in [0,1]
 * to two partitions of a graph, by considering pairs of nodes.
 * This implementation approximates the index by sampling.
 */
class SampledNodeStructuralRandMeasure: public NetworKit::DissimilarityMeasure {

public:

	/**
	 * Constructs the SampledNodeStructuralRandMeasure. A maximum of @a maxSamples samples are drawn.
	 *
	 * @param maxSamples The amount of samples to draw.
	 */
	SampledNodeStructuralRandMeasure(count maxSamples);

	virtual double getDissimilarity(const Graph& G, const Partition& first, const Partition& second);

protected:

	count maxSamples;

};

} /* namespace NetworKit */
#endif /* SAMPLEDRANDMEASURE_H_ */
