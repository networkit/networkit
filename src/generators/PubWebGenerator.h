/*
 * PubWebGenerator.h
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#ifndef PUBWEBGENERATOR_H_
#define PUBWEBGENERATOR_H_

#include "StaticGraphGenerator.h"

#define _USE_MATH_DEFINES // to be able to use M_PI
#include <cmath>



namespace NetworKit {

const int MAX_RAND_VAL = 1000;
const float MAX_DENSE_AREA_RADIUS = 0.35f;
const float MIN_MAX_DENSE_AREA_FACTOR = 5.0f;


struct circle {
	float x;
	float y;
	float rad;
};

/**
 * Generates a static graph that resembles an assumed geometric distribution of nodes in
 * a P2P network. The basic structure is to distribute points randomly in the unit torus
 * and to connect vertices close to each other (at most @a neighRad distance and none of
 * them already has @a maxNeigh neighbors). The distribution is chosen to get some areas with
 * high density and others with low density. There are @a numDenseAreas dense areas, which can
 * overlap. Each area is circular, has a certain position and radius and number of points.
 * These values are strored in @a denseAreaXYR and @a numPerArea, respectively.
 *
 * Used and described in more detail in J. Gehweiler, H. Meyerhenke: A Distributed
 * Diffusive Heuristic for Clustering a Virtual P2P Supercomputer. In Proc. 7th High-Performance
 * Grid Computing Workshop (HPGC'10), in conjunction with 24th IEEE Internatl. Parallel and
 * Distributed Processing Symposium (IPDPS'10), IEEE, 2010.
 *
 * Reasonable parameters for constructor:
 * - numNodes: up to a few thousand (possibly more if visualization is not desired and quadratic
 *   time complexity has been resolved)
 * - numberOfDenseAreas: [10, 50]
 * - neighborhoodRadius: [0.1, 0.35]
 * - maxNumberOfNeighbors: [4, 40]
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
	PubWebGenerator() {} // nullary constructor needed for Python Shell - do not use this to construct instance

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
