/*
 * ErdosRenyiGenerator.h
 *
 *  Created on: 21.01.2014
 *      Author: Henning
 */

#ifndef ERDOSRENYIGENERATOR_H_
#define ERDOSRENYIGENERATOR_H_

#include "StaticGraphGenerator.h"
#include <cmath>

namespace NetworKit {

/**
 * Creates G(n, p) graphs.
 */
class ErdosRenyiGenerator: public StaticGraphGenerator {
protected:
	count n;
	double p;

public:
	ErdosRenyiGenerator(count nNodes, double prob);
	virtual ~ErdosRenyiGenerator();

	virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* ERDOSRENYIGENERATOR_H_ */
