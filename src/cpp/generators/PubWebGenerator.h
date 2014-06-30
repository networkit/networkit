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
const float MAX_DENSE_AREA_RADIUS = 0.2f;
const float MIN_MAX_DENSE_AREA_FACTOR = 5.0f;
const edgeweight BASE_WEIGHT = 0.01f;

typedef float distance; // TODO: use double, not float
typedef std::pair<node, node> edge;


struct circle {
	float x;
	float y;
	float rad;
};

/**
 * @ingroup generators
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
 * - numNodes: from several hundred up to a few thousand
 *   (possibly more if visualization is not desired and quadratic time
 *   complexity has been resolved)
 * - numberOfDenseAreas: depending on number of nodes, e.g. [8, 50]
 * - neighborhoodRadius: the higher, the better the connectivity [0.1, 0.35]
 * - maxNumberOfNeighbors: maximum degree, a higher value corresponds to better connectivity [4, 40]
 */
class PubWebGenerator: public NetworKit::StaticGraphGenerator {

	friend class DynamicPubWebGenerator;

protected:
	count n; //!< number of nodes
	count numDenseAreas; //!< number of areas with more nodes (denser)
	float neighRad; //!< neighborhood radius
	count maxNeigh; //!< maximum number of neighbors
	std::vector<circle> denseAreaXYR; //!< position of each circular dense area
	std::vector<count> numPerArea; //!< number of points in each circular area

	void determineNeighbors(Graph& g);
	void moveNodeIntoUnitSquare(float& x, float& y);
	float squaredDistanceInUnitTorus(float x1, float y1, float x2, float y2);
	void chooseDenseAreaSizes();
	void fillDenseAreas(Graph& g);
	void spreadRemainingNodes(Graph& g);
	void chooseClusterSizes();
	void addNodesToArea(index area, count num, Graph& g);
	bool isValidEdge(Graph& g, node u, node v, edgeweight& ew);

public:
	PubWebGenerator() {} // nullary constructor needed for Python Shell - do not use this to construct instance

	PubWebGenerator(count numNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors);

	virtual Graph generate();

protected:

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
