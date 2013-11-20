/*
 * PubWebGenerator.h
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#ifndef PUBWEBGENERATOR_H_
#define PUBWEBGENERATOR_H_

#include "StaticGraphGenerator.h"
#include "../auxiliary/RandomProbability.h"
#include "../auxiliary/RandomInteger.h"
#include <cmath>

namespace NetworKit {

struct circle {
	float x;
	float y;
	float rad;
};

/**
 * TODO: class description
 */
class PubWebGenerator: public NetworKit::StaticGraphGenerator {

	friend class DynamicPubWebGenerator;

protected:
	count n; //!< number of nodes
	count numDenseAreas; //!< number of areas with more nodes (denser)
	float neighRad; //!< neighborhood radius
	count maxNeigh; //!< maximum number of neighbors
	std::vector<circle> denseAreaXYR;
	std::vector<count> numPerArea;

	void determineNeighbors(Graph& g);
	void determineNeighborsOf(Graph& g, node u);
	void moveNodeIntoUnitSquare(float& x, float& y);
	float squaredDistanceInUnitTorus(float x1, float y1, float x2, float y2);
	void chooseDenseAreaSizes();
	void fillDenseAreas(Graph& g);
	void spreadRemainingNodes(Graph& g);
	void chooseClusterSizes();
	void addNodesToArea(index area, count num, Graph& g);
	bool isValidEdge(Graph& g, node u, node v);

public:
	PubWebGenerator(); // nullary constructor needed for Python Shell - do not use this to construct instance

	// TODO: provide reasonable default parameters in documentation or as default arguments to constructor
	PubWebGenerator(count numNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors);
	virtual ~PubWebGenerator();

	virtual Graph generate();

	/**
	 * Adds nodes randomly, distribution respects original one.
	 */
	void addNode(Graph& g);

	/**
	 * Removes random node, uniform distribution.
	 */
	void removeRandomNode(Graph& g);
};

} /* namespace NetworKit */
#endif /* PUBWEBGENERATOR_H_ */
