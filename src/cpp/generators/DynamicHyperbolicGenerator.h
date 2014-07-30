/*
 * DynamicHyperbolicGenerator.h
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#ifndef DYNAMICHYPERBOLICGENERATOR_H_
#define DYNAMICHYPERBOLICGENERATOR_H_

#include "DynamicGraphGenerator.h"
#include "Quadtree/Quadtree.h"

namespace NetworKit {

class DynamicHyperbolicGenerator: public NetworKit::DynamicGraphGenerator  {
public:
	DynamicHyperbolicGenerator(count n, double initialFactor = 1, double alpha = 1, double stretch = 1, double moveEachStep = 0, double factorgrowth = 0, double moveDistance = 0);
	DynamicHyperbolicGenerator(std::vector<double> &angles, std::vector<double> &radii, double R, double initialFactor = 1, double moveEachStep = 0, double factorgrowth = 0, double moveDistance = 0);
	DynamicHyperbolicGenerator();
	virtual ~DynamicHyperbolicGenerator();
	std::vector<GraphEvent> generate(count nSteps) override;
	virtual void initializeGraph();

private:
	double factorgrowth;
	double moveDistance;
	count nodes;
	double moveEachStep;
	double stretch;
	double currentfactor;
	double alpha;
	Quadtree<index> quad;
	vector<double> angles;
	vector<double> radii;
	bool initialized;
};

} /* namespace NetworKit */
#endif /* DYNAMICHYPERBOLICGENERATOR_H_ */

