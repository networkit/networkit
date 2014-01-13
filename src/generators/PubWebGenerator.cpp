/*
 * PubWebGenerator.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include "../auxiliary/Random.h"

#include "PubWebGenerator.h"

namespace NetworKit {


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

}

PubWebGenerator::~PubWebGenerator() {

}

bool PubWebGenerator::isValidEdge(Graph& g, node u, node v) {
	auto isValid([&](node u, node v, float squaredDistance) {
		return ((squaredDistance <= neighRad * neighRad)
				&& (g.degree(u) < maxNeigh)
				&& (g.degree(v) < maxNeigh));
	});

	float x1 = g.getCoordinate(u, 0);
	float y1 = g.getCoordinate(u, 1);
	float x2 = g.getCoordinate(v, 0);
	float y2 = g.getCoordinate(v, 1);
	float sqrDist = squaredDistanceInUnitTorus(x1, y1, x2, y2);

	return isValid(u, v, sqrDist);
}

void PubWebGenerator::determineNeighborsOf(Graph& g, node u) {
	g.forNodes([&](node v) {
		if (isValidEdge(g, u, v)) {
			g.addEdge(u, v);
		}
	});
}

void PubWebGenerator::determineNeighbors(Graph& g) {
	g.forNodePairs([&](node u, node v) { // TODO: improve quadratic loop!
		if (isValidEdge(g, u, v)) {
			g.addEdge(u, v);
		}
	});
}

void PubWebGenerator::addNodesToArea(index area, count num, Graph& g) {

	for (index j = 0; j < num; ++j) {
		// compute random angle between [0, 2pi) and distance between [0, width/2]
		float angle = Aux::Random::real() * 2.0 * M_PI;
		float dist = Aux::Random::real() * denseAreaXYR[area].rad;

		// compute coordinates and adjust them
		float x = denseAreaXYR[area].x + cosf(angle) * dist;
		float y = denseAreaXYR[area].y + sinf(angle) * dist;
		moveNodeIntoUnitSquare(x, y);

		// create vertex with these coordinates
		g.addNode(x, y);
	}
}

void PubWebGenerator::fillDenseAreas(Graph& g) {

	for (index area = 0; area < numDenseAreas; ++area) {
		// choose center randomly, ensure complete cluster is within (0,1) without modifications
		denseAreaXYR[area].x = Aux::Random::real();
		denseAreaXYR[area].y = Aux::Random::real();
		addNodesToArea(area, numPerArea[area], g);
	}
}

void PubWebGenerator::chooseDenseAreaSizes() {
	denseAreaXYR.reserve(numDenseAreas);

	for (index area = 0; area < numDenseAreas; ++area) {
		// anti-quadratic probability distribution
		float f = Aux::Random::real() * MIN_MAX_DENSE_AREA_FACTOR + 1.0f;
		denseAreaXYR[area].rad = (MAX_DENSE_AREA_RADIUS * f * f)
				/ (MIN_MAX_DENSE_AREA_FACTOR * MIN_MAX_DENSE_AREA_FACTOR);
	}
}

// randomly spread remaining vertices over whole area
void PubWebGenerator::spreadRemainingNodes(Graph& g) {

	while (g.numberOfNodes() < n) {
		float x = Aux::Random::real();
		float y = Aux::Random::real();
		g.addNode(x, y);
	}
}

// compute number of nodes per cluster, each cluster has approx. same density
void PubWebGenerator::chooseClusterSizes() {
	float f = 0.0;
	for (index i = 0; i < numDenseAreas; ++i) {
		f += pow(denseAreaXYR[i].rad, 1.5);
	}
	f = ((float) n * ((float) numDenseAreas / ((float) numDenseAreas + 2.0f)))
			/ f;
	// TODO: better formula?

	numPerArea.reserve(numDenseAreas);
	for (index i = 0; i < numDenseAreas; ++i) {
		numPerArea[i] = roundf(f * pow(denseAreaXYR[i].rad, 1.5));
	}
}

Graph PubWebGenerator::generate() {
	// init
	Graph g(0);

	// add vertices according to PubWeb distribution
	chooseDenseAreaSizes();
	chooseClusterSizes();
	fillDenseAreas(g);
	spreadRemainingNodes(g);
	determineNeighbors(g);

	return g;
}


// TODO: NOT tested!
void PubWebGenerator::removeRandomNode(Graph& g) {
	node u = Aux::Random::integer(n - 1);
	g.removeNode(u);
}

} /* namespace NetworKit */

