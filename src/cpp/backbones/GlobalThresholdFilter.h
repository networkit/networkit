/*
 * GlobalThresholdFilter.h
 *
 *  Created on: 23.07.2014
 *      Author: Gerd Lindner
 */

#ifndef GLOBALTHRESHOLDFILTER_H_
#define GLOBALTHRESHOLDFILTER_H_

#include "AttributeGenerator.h"
#include "BackboneCalculator.h"

namespace NetworKit {

/** 
 * Calculates a backbone by applying a global threshold to an edge attribute.
 */
class GlobalThresholdFilter : public BackboneCalculator {

public:

	/**
	 * Creates a new instance of a global threshold filter.
	 * @param threshold		the threshold
	 * @param above			if set to true, edge attribute needs to be above or equal to the threshold.
	 * 						If set to false, edge attribute needs to be below or equal to the threshold.
	 */
	GlobalThresholdFilter(double threshold, bool above); //TODO: better name for parameter?

	Graph calculate(const Graph& graph, const EdgeAttribute& attribute);

private:
	double threshold;
	bool above;

};

}
/* namespace NetworKit */
#endif /* GLOBALTHRESHOLDFILTER_H_ */
