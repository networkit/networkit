/*
 * DynamicPubWebGenerator.h
 *
 *  Created on: 15.01.2014
 *      Author: Henning
 */

#ifndef DYNAMICPUBWEBGENERATOR_H_
#define DYNAMICPUBWEBGENERATOR_H_

#include "DynamicGraphGenerator.h"
#include "PubWebGenerator.h"
#include "../dynamics/GraphEvent.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

class DynamicPubWebGenerator: public NetworKit::DynamicGraphGenerator
{

protected:
	PubWebGenerator initGen; // multiple inheritance did not work with different generate functions

public:
	DynamicPubWebGenerator(count numNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors);

	virtual ~DynamicPubWebGenerator();

	/**
	 * Generate event stream.
	 *
	 * @param[in]	nSteps	number of time steps in the event stream
	 */
	virtual std::vector<GraphEvent> generate(count nSteps);
};

} /* namespace NetworKit */
#endif /* DYNAMICPUBWEBGENERATOR_H_ */
