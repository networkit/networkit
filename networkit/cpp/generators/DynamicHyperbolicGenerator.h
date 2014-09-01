/*
 * DynamicHyperbolicGenerator.h
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#ifndef DYNAMICHYPERBOLICGENERATOR_H_
#define DYNAMICHYPERBOLICGENERATOR_H_

#include <map>

#include "DynamicGraphGenerator.h"
#include "Quadtree/Quadtree.h"


namespace NetworKit {

class DynamicHyperbolicGenerator: public NetworKit::DynamicGraphGenerator  {
	friend class GeneratorsGTest;
public:
	DynamicHyperbolicGenerator(count n, double initialFactor = 1, double alpha = 1, double stretch = 1, double moveEachStep = 0, double factorgrowth = 0, double moveDistance = 0);
	DynamicHyperbolicGenerator(std::vector<double> &angles, std::vector<double> &radii, double stretch, double initialFactor = 1, double moveEachStep = 0, double factorgrowth = 0, double moveDistance = 0);
	DynamicHyperbolicGenerator();
	virtual ~DynamicHyperbolicGenerator();
	std::vector<GraphEvent> generate(count nSteps) override;
	void initializeQuadTree();
	void initializeMovement();
	Graph getGraph();



	/** only to be used by cython - this eliminates an unnecessary copy */
	Graph* _getGraph() {
		return new Graph{std::move(getGraph())};
	};

	std::vector<Point<float> > getCoordinates() const;
	std::vector<Point<float> > getHyperbolicCoordinates() const;

private:
	count quadTreeHeight() {
		return quad.height();
	}
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
	vector<double> angularMovement;
	vector<double> radialMovement;
	bool initialized;
};

} /* namespace NetworKit */
#endif /* DYNAMICHYPERBOLICGENERATOR_H_ */

