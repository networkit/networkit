/*
 * SimmelianMultiscaleBackbone.h
 *
 *  Created on: 16.07.2014
 *      Author: Gerd Lindner
 */

#ifndef SIMMELIANMULTISCALEBACKBONE_H_
#define SIMMELIANMULTISCALEBACKBONE_H_

#include "BackboneCalculator.h"
#include "gtest/gtest_prod.h"
#include <set>

namespace NetworKit {


/** 
 * Combines triangle counting and multiscale backbone filtering.
 */
class SimmelianMultiscaleBackbone : public BackboneCalculator {

public:

	/**
	 * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
	 * @param threshold 	the filtering value in (0,1) that is used in multiscale filtering.
	 */
	SimmelianMultiscaleBackbone(double threshold);

	Graph calculate(const Graph& graph);

private:
	//Calculation parameters
	double threshold;
};

}
/* namespace NetworKit */
#endif /* SIMMELIANMULTISCALEBACKBONE_H_ */
