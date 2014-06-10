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
	/**
	 * Calculates the backbone graph for the given input graph.
	 */
	virtual Graph calculate() = 0;

	/** only to be used by cython - this eliminates an unnecessary copy */
	Graph* _calculate() {
		return new Graph{std::move(calculate())};
	};
};

} /* namespace NetworKit */
#endif /* BACKBONECALCULATOR_H_ */
