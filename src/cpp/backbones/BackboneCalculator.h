/*
 * BackboneCalculator.h
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#ifndef BACKBONECALCULATOR_H_
#define BACKBONECALCULATOR_H_

#include "../graph/Graph.h"

namespace NetworKit {

/** 
 * Abstract base class for Backbone Calculators.
 */
class BackboneCalculator {

public:
	BackboneCalculator();
	virtual ~BackboneCalculator();

	/**
	 * Calculates the backbone graph for the given input graph.
	 */
	virtual Graph calculate(const Graph& g) = 0;

	virtual count test() = 0;
};

} /* namespace NetworKit */
#endif /* BACKBONECALCULATOR_H_ */
