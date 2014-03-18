/*
 * RmatGenerator.h
 *
 *  Created on: 18.03.2014
 *      Author: Henning
 */

#ifndef RMATGENERATOR_H_
#define RMATGENERATOR_H_

#include "StaticGraphGenerator.h"
#include "../graph/Graph.h"

namespace NetworKit {

class RmatGenerator: public NetworKit::StaticGraphGenerator {
protected:
	count scale; ///< n = 2^scale
	count edgeFactor;
	double a, b, c, d; ///< probabilities
	double defaultEdgeWeight;

public:
	/**
	 * DO NOT use to create instances. For Python interface only.
	 */
	RmatGenerator();

	/**
	 * @param[in] scale Number of nodes = 2^scale
	 * @param[in] edgeFactor Number of edges = number of nodes * edgeFactor
	 * @param[in] a Probability for quadrant upper left
	 * @param[in] b Probability for quadrant upper right
	 * @param[in] c Probability for quadrant lower left
	 * @param[in] d Probability for quadrant lower right
	 */
	RmatGenerator(count scale, count edgeFactor, double a, double b, double c, double d);
	virtual ~RmatGenerator();

	Graph generate() override;
};

} /* namespace NetworKit */
#endif /* RMATGENERATOR_H_ */
