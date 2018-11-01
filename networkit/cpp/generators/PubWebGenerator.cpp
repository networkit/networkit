/*
 * PubWebGenerator.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include <set>
#include <queue>

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

bool PubWebGenerator::isValidEdge(Graph& g, node u, node v, edgeweight& weight) {

	auto isValid([&](node u, node v, float squaredDistance) {
		return ((squaredDistance <= neighRad * neighRad)
				&& (g.degree(u) < maxNeigh)
				&& (g.degree(v) < maxNeigh));
	});

	Point<float> pu = g.getCoordinate(u);
	Point<float> pv = g.getCoordinate(v);
	float& x1 = pu[0];
	float& y1 = pu[1];
	float& x2 = pv[0];
	float& y2 = pv[1];
	float sqrDist = squaredDistanceInUnitTorus(x1, y1, x2, y2);

	weight = BASE_WEIGHT /  sqrt(sqrDist);

	return isValid(u, v, sqrDist);
}


// TODO: use ANN or similar library with appropriate space-partitioning data structure to
//       get rid of quadratic time complexity
void PubWebGenerator::determineNeighbors(Graph& g) {

	float sqrNeighRad = neighRad * neighRad;
	std::set<edge> eligibleEdges;

	auto isInRange([&](float squaredDistance) {
		return (squaredDistance <= sqrNeighRad);
	});

	g.forNodes([&](node u) {
		std::priority_queue<std::pair<distance, edge> > pq;
		Point<float> p1 = g.getCoordinate(u);
		float& x1 = p1[0];
		float& y1 = p1[1];

		// fill PQ with neighbors in range
		g.forNodes([&](node v) {
			Point<float> p2 = g.getCoordinate(v);
			float& x2 = p2[0];
			float& y2 = p2[1];
			float sqrDist = squaredDistanceInUnitTorus(x1, y1, x2, y2);

			if (isInRange(sqrDist)) {
				edge e = std::make_pair(std::min(u, v), std::max(u, v));
				pq.push(std::make_pair(-sqrDist, e));
			}
		});

		// mark maxNeigh nearest neighbors as eligible or insert them into graph g
		count end = std::min(maxNeigh, (count) pq.size());
		for (index i = 0; i < end; ++i) {
			std::pair<distance, edge> currentBest = pq.top();
			pq.pop();

			if (eligibleEdges.count(currentBest.second) > 0) {
				// edge is already marked => insert it
				edgeweight ew = BASE_WEIGHT / -currentBest.first;
				g.addEdge(currentBest.second.first, currentBest.second.second, ew);
//				TRACE("add edge " , currentBest.second.first , "/" , currentBest.second.second);
			}
			else {
				// mark edge as eligible
				eligibleEdges.insert(currentBest.second);
			}
		}
	});
}

void PubWebGenerator::addNodesToArea(index area, count num, Graph& g) {

	for (index j = 0; j < num; ++j) {
		// compute random angle between [0, 2pi) and distance between [0, width/2]
		float angle = Aux::Random::real() * 2.0 * PI;
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
	denseAreaXYR.resize(numDenseAreas);

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

	numPerArea.resize(numDenseAreas);
	for (index i = 0; i < numDenseAreas; ++i) {
		numPerArea[i] = roundf(f * pow(denseAreaXYR[i].rad, 1.5));
	}
}

Graph PubWebGenerator::generate() {
	// init
	Graph G(0, true);

	// add vertices according to PubWeb distribution
	chooseDenseAreaSizes();
	chooseClusterSizes();
	fillDenseAreas(G);
	spreadRemainingNodes(G);
	determineNeighbors(G);

	G.shrinkToFit();
	return G;
}


// TODO: NOT tested!
void PubWebGenerator::removeRandomNode(Graph& g) {
	node u = Aux::Random::integer(n - 1);
	g.removeNode(u);
}

} /* namespace NetworKit */

