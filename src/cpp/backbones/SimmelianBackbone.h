/*
 * SimmelianBackbone.h
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#ifndef SIMMELIANBACKBONE_H_
#define SIMMELIANBACKBONE_H_

#include "BackboneCalculator.h"

namespace NetworKit {

/** 
 * Calculates the simmelian backbone for a given input graph.
 */
class SimmelianBackbone : public BackboneCalculator {

public:
	SimmelianBackbone();
	~SimmelianBackbone();

	/**
	 * Calculates the simmelian backbone for the given graph.
	 */
	virtual Graph calculate(const Graph& g);

	virtual count test();
};

} /* namespace NetworKit */
#endif /* SIMMELIANBACKBONE_H_ */
