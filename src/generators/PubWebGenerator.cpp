/*
 * PubWebGenerator.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include "PubWebGenerator.h"

namespace NetworKit {

// FIXME: change to constants
#define THRSH_CENTER_ITER 10
#define ONE_DIM_SCALE 0.0f
#define MAX_RAND_VAL 1000
#define BENEFIT_INC 0.25f
#define PI 3.141f
#define MAX_DENSE_AREA_RADIUS 0.25f
#define MIN_MAX_DENSE_AREA_FACTOR 5 // ###
#define PART_CHANGE_FACTOR 0.0001f
#define INIT_LOAD 100.0f

void PubWebGenerator::moveNodeIntoUnitSquare(float& x, float& y) {
	auto move([&](float& z) {
		if (z > 1.0f) {
			z -= 1.0f;
		} else if (z < 0.0f) {
			z += 1.0f;
		}
	});

	move(x);
	move(y);
}

float PubWebGenerator::squaredDistanceInUnitTorus(float x1, float y1, float x2,
		float y2) {
	auto adjustForUnitTorus([&](float& z) {
		if (z > 0.5) {
			z = 1.0 - z;
		}
		else if (z < -0.5) {
			z = z + 1.0;
		}
	});

	float distx = x1 - x2;
	float disty = y1 - y2;
	adjustForUnitTorus(distx);
	adjustForUnitTorus(disty);

	return (distx * distx + disty * disty);
}

PubWebGenerator::PubWebGenerator(count numNodes, count numberOfDenseAreas,
		float neighborhoodRadius, count maxNumberOfNeighbors) :
		n(numNodes), numDenseAreas(numberOfDenseAreas), neighRad(
				neighborhoodRadius), maxNeigh(maxNumberOfNeighbors) {
	// TODO Auto-generated constructor stub

}

PubWebGenerator::~PubWebGenerator() {
	// TODO Auto-generated destructor stub
}

void PubWebGenerator::determineNeighbors(Graph& g) {
	auto isValidEdge([&](node u, node v, float squaredDistance) {
		return ((squaredDistance <= neighRad * neighRad)
				&& (g.degree(u) <= maxNeigh)
				&& (g.degree(v) <= maxNeigh));
	});

	g.forNodePairs([&](node u, node v) { // FIXME: quadratic loop!
				float sqrDist = squaredDistanceInUnitTorus(g.getCoordinate(u, 0), g.getCoordinate(u, 1), g.getCoordinate(v, 0), g.getCoordinate(v, 1));
				if (isValidEdge(u, v, sqrDist)) {
					g.addEdge(u, v);
				}
			});
}

std::vector<float> PubWebGenerator::fillDenseAreas(Graph& g,
		const std::vector<count>& numPerArea,
		std::vector<circle>& denseAreaXYR) {
	std::vector<float> coordinates;
	Aux::RandomProbability randGen;

	for (index area = 0; area < numDenseAreas; ++area) {
		// choose center randomly, ensure complete cluster is within (0,1) without modifications
		denseAreaXYR[area].x = randGen.randomFloat();
		denseAreaXYR[area].y = randGen.randomFloat();

		for (index j = 0; j < numPerArea[area]; ++j) {
			// compute random angle between [0, 2pi) and distance between [0, width/2]
			float angle = randGen.randomFloat() * 2.0 * PI;
			float dist = randGen.randomFloat() * denseAreaXYR[area].rad;

			// compute coordinates and adjust them
			float x = denseAreaXYR[area].x + cosf(angle) * dist;
			float y = denseAreaXYR[area].y + sinf(angle) * dist;
			moveNodeIntoUnitSquare(x, y);

			// create vertex with this coordinate
			g.addNode();
			coordinates.push_back(x);
			coordinates.push_back(y);
		}
	}

	return coordinates;
}

std::vector<circle> PubWebGenerator::chooseDenseAreaSizes() {
	std::vector<circle> denseAreaXYR(numDenseAreas);
	Aux::RandomProbability randGen;

	for (index area = 0; area < numDenseAreas; ++area) {
		// anti-quadratic probability distribution
		float f = randGen.randomFloat() * MIN_MAX_DENSE_AREA_FACTOR + 1.0f;
		denseAreaXYR[area].rad = (MAX_DENSE_AREA_RADIUS * f * f)
				/ (MIN_MAX_DENSE_AREA_FACTOR * MIN_MAX_DENSE_AREA_FACTOR);
	}

	return denseAreaXYR;
}

// randomly spread remaining vertices over whole area
void PubWebGenerator::spreadRemainingNodes(Graph& g, std::vector<float>& coordinates) {
	Aux::RandomProbability randGen;

	while (coordinates.size() < 2 * n) {
		g.addNode();
		coordinates.push_back(randGen.randomFloat());
		coordinates.push_back(randGen.randomFloat());
	}
}

// compute number of nodes per cluster, each cluster has approx. same density
std::vector<count> PubWebGenerator::chooseClusterSizes(std::vector<circle>& denseAreaXYR) {
	float f = 0.0;
	for (index i = 0; i < numDenseAreas; ++i) {
		f += pow(denseAreaXYR[i].rad, 1.5);
	}
	f = ((float) n * ((float) numDenseAreas / ((float) numDenseAreas + 2.0f))) / f;
	// TODO: better formula?

	std::vector<count> numPerArea(numDenseAreas);
	for (index i = 0; i < numDenseAreas; ++i) {
		numPerArea[i] = roundf(f * pow(denseAreaXYR[i].rad, 1.5));
	}

	return numPerArea;
}

Graph PubWebGenerator::generate() {
	// init
	Graph g(0);
	count dims = 2;
	std::vector<float> coordinates;
	Aux::RandomProbability randGen;

	// choose area sizes
	std::vector<circle> denseAreaXYR = chooseDenseAreaSizes();

	// compute number of nodes per cluster, each cluster has approx. same density
	std::vector<count> numPerArea = chooseClusterSizes(denseAreaXYR);

	// fill dense areas
	coordinates = fillDenseAreas(g, numPerArea, denseAreaXYR);

	// randomly spread remaining vertices over whole area
	spreadRemainingNodes(g, coordinates);

	// insert coordinates into graph
	g.initCoordinates(dims);
	for (index v = 0; v < n; ++v) {
		g.setCoordinate(v, 0, coordinates[v * dims]);
		g.setCoordinate(v, 1, coordinates[v * dims + 1]);
	}

	// determine neighbors before adjusting the coordinates
	determineNeighbors(g);

	// adjust coordinates for postscript output
	g.forNodes([&](node u) {
		g.setCoordinate(u, 0, 1000.0 * g.getCoordinate(u, 0));
		g.setCoordinate(u, 1, 1000.0 * g.getCoordinate(u, 1));
	});

	return g;
}

} /* namespace NetworKit */

