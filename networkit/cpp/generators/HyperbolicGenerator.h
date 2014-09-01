/*
 * HyperbolicGenerator.h
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef HYPERBOLICGENERATOR_H_
#define HYPERBOLICGENERATOR_H_

#include <vector>
#include <map>
#include "../geometric/HyperbolicSpace.h"
#include "StaticGraphGenerator.h"

namespace NetworKit {

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based

class HyperbolicGenerator: public NetworKit::StaticGraphGenerator {
public:
	HyperbolicGenerator();
	HyperbolicGenerator(count n, double factor = 1, double alpha = 1, double stretchradius = 1);
	HyperbolicGenerator(count n, count m);
	virtual ~HyperbolicGenerator();
	static double expectedNumberOfEdges(count n, double distanceFactor, double stretch);
	static Graph generate(vector<double> * angles, vector<double> * radii, double R, double thresholdDistance);
	static std::map<index, Point<float> > getCoordinates(vector<double> &angles, vector<double> &radii);
	Graph generate(count n, double distanceFactor=1, double alpha=1, double stretchradius = 1);
	Graph generate();

private:
	count nodeCount;
	double stretch;
	double factor;
	double alpha;

};
}
#endif /* HYPERBOLICGENERATOR_H_ */
