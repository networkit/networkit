/*
 * HyperbolicGenerator.h
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
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
	Graph generate(vector<double> * angles, vector<double> * radii, double R, double thresholdDistance);
	Graph generate(count n, double stretchradius = 1);//TODO: add converter to NetworKit graph
	Graph generate();

	count nodeCount;
	double stretch;
};
}
#endif /* HYPERBOLICGENERATOR_H_ */
