/*
 * HyperbolicGenerator.h
 *
 *  Created on: 20.05.2014
 *      Author: moritz
 */

#ifndef HYPERBOLICGENERATOR_H_
#define HYPERBOLICGENERATOR_H_

#include <vector>
#include <random>
#include "HyperbolicSpace.h"
#include "StaticGraphGenerator.h"

namespace NetworKit {

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based

class HyperbolicGenerator: public NetworKit::StaticGraphGenerator {
public:
	HyperbolicGenerator();
	HyperbolicGenerator(count n, double stretchradius);
	virtual ~HyperbolicGenerator();
	Graph generate(count n, double stretchradius = 1);//TODO: add converter to NetworKit graph
	Graph generate();

	count nodeCount;
	double stretch;
};
}
#endif /* HYPERBOLICGENERATOR_H_ */
