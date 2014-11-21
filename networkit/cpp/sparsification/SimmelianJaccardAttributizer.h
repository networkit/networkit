/*
 * SimmelianJaccardAttributizer.h
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#ifndef SIMMELIANJACCARDATTRIBUTIZER_H_
#define SIMMELIANJACCARDATTRIBUTIZER_H_

#include "SimmelianAttributizer.h"

namespace NetworKit {

/**
 * * Calculates the Simmelian backbone (non-parametric aka jaccard variant) for a given input graph.
 */
class SimmelianJaccardAttributizer : public SimmelianAttributizer {

public:

	/**
	 * Creates a new instance of the non-parametric (jaccard) variant of the Simmelian Backbone calculator.
	 */
	SimmelianJaccardAttributizer(const Graph& graph, const std::vector<int>& triangles);

	virtual std::vector<double> getAttribute() override;

};

}
/* namespace NetworKit */

#endif /* SIMMELIANJACCARDATTRIBUTIZER_H_ */
