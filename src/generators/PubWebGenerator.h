/*
 * PubWebGenerator.h
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#ifndef PUBWEBGENERATOR_H_
#define PUBWEBGENERATOR_H_

#include "StaticGraphGenerator.h"
#include "../aux/RandomProbability.h"
#include <cmath>

namespace NetworKit {

struct circle {
	float x;
	float y;
	float rad;
};

class PubWebGenerator: public NetworKit::StaticGraphGenerator {
protected:
	count n; //!< number of nodes
	count numDenseAreas; //!< number of areas with more nodes (denser)
	float neighRad; //!< neighborhood radius
	count maxNeigh; //!< maximum number of neighbors
	std::vector<circle> denseAreaXYR;

	void determineNeighbors(Graph& g);
	void moveNodeIntoUnitSquare(float& x, float& y);
	float squaredDistanceInUnitTorus(float x1, float y1, float x2, float y2);
	void chooseDenseAreaSizes();
	std::vector<float> fillDenseAreas(Graph& g,
			const std::vector<count>& numPerArea);
	void spreadRemainingNodes(Graph& g, std::vector<float>& coordinates);
	std::vector<count> chooseClusterSizes();

public:
	PubWebGenerator(count numNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors);
	virtual ~PubWebGenerator();

	virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* PUBWEBGENERATOR_H_ */
