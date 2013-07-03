/*
 * SampledRandMeasure.h
 *
 *  Created on: 01.07.2013
 *      Author: cls
 */

#ifndef SAMPLEDRANDMEASURE_H_
#define SAMPLEDRANDMEASURE_H_

#include "DissimilarityMeasure.h"

#include "../auxiliary/RandomInteger.h"

namespace NetworKit {

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
