/*
 * DynamicPubWebGenerator.cpp
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */


#include "DynamicPubWebGenerator.h"

namespace NetworKit {



DynamicPubWebGenerator::DynamicPubWebGenerator(GraphEventProxy& proxy,
		count numInitialNodes, count numberOfDenseAreas, float neighborhoodRadius,
		count maxNumberOfNeighbors) :
		DynamicGraphGenerator(proxy)
//,
//		staticGen(numInitialNodes, numberOfDenseAreas, neighborhoodRadius,
//				maxNumberOfNeighbors)
{

}

DynamicPubWebGenerator::~DynamicPubWebGenerator() {
	// TODO Auto-generated destructor stub
}


/**
 * Assumption: Static PubWeb graph already exists.
 */
void DynamicPubWebGenerator::initializeGraph() {
	// add vertices according to PubWeb distribution
	// TODO: disable comment
//	this->chooseDenseAreaSizes();
//	this->chooseClusterSizes();
//	this->fillDenseAreas(*G);
//	this->spreadRemainingNodes(*G);
//	this->determineNeighbors(*G);
}


#if 0
void Move(int V, float vmin, float vmax, int stop, float * xy, float * way, float * vel, int * vstop)
{
  int i;
  float distx, disty, dist;

  for (i = 0; i < V; i++) {
    if (vstop[i] > 0)
      // Pause
      vstop[i]--;
    else {
      distx = way[i + i] - xy[i + i];
      disty = way[i + i + 1] - xy[i + i + 1];

      // Wrap around
      AdjustWrapAround(distx, disty);
      dist = sqrt(distx * distx + disty * disty);

      if (dist < vel[i]) {
	vstop[i] = rand() % stop;
	xy[i + i] = way[i + i];
	xy[i + i + 1] = way[i + i + 1];
	// Neues Ziel
	way[i + i] = (float)rand() / (float)RAND_MAX;
	way[i + i + 1] = (float)rand() / (float)RAND_MAX;
#ifdef ONE_DIM
	way[i + i + 1] *= ONE_DIM_SCALE;
#endif
	vel[i] = vmin + ((float)rand() / (float)RAND_MAX) * (vmax - vmin);
      }
      else { // berechne Zwischenhalt
	xy[i + i] += distx * vel[i] / dist;

	if (xy[i + i] > 1.0)
	  xy[i + i] -= 1.0;
	if (xy[i + i] < 0.0)
	  xy[i + i] += 1.0;

	xy[i + i + 1] += disty * vel[i] / dist;

	if (xy[i + i + 1] > 1.0)
	  xy[i + i + 1] -= 1.0;
	if (xy[i + i + 1] < 0.0)
	  xy[i + i + 1] += 1.0;
      }
    }
  }
}
#endif


void DynamicPubWebGenerator::moveNodesRandomly() {
#if 0
	Aux::RandomProbability randGen;

	this->G->forNodes([&](node u) {
		// current position
		float x = this->G->getCoordinate(u, 0);
		float y = this->G->getCoordinate(u, 1);

		// compute random direction
		float angle = randGen.randomFloat() * M_PI_2;

		// compute random distance
		float dist = randGen.randomFloat();

		// compute new positions
		x += cosf(angle) * dist;
		y += sinf(angle) * dist;

		// adjust for wraparound
		this->moveNodeIntoUnitSquare(x, y);

		// move node
		this->Gproxy->G->setCoordinate(u, 0, x);
		this->Gproxy->G->setCoordinate(u, 1, y);
	});
#endif
}


void DynamicPubWebGenerator::generate() {
#if 0
	// nodes move to new position
	moveNodesRandomly();

	// decide if new edges have to be inserted
	// FIXME: inefficient!
	G->forNodePairs([&](node u, node v) {
		bool valid = this->isValidEdge(*G, u, v);
		bool exists = G->hasEdge(u, v);

		if (valid && ! exists) {
			// insert
			this->Gproxy->addEdge(u, v);
		}
		if (! valid && exists) {
			// remove
			this->Gproxy->removeEdge(u, v);
		}
	});

	this->Gproxy->timeStep(); // trigger a time step
#endif
}


#if 0


/** FIXME: the following is code duplicated from the static generator
 *  Reason: Need to use the proxy
 **/

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

void DynamicPubWebGenerator::moveNodeIntoUnitSquare(float& x, float& y) {
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

float DynamicPubWebGenerator::squaredDistanceInUnitTorus(float x1, float y1, float x2,
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

bool DynamicPubWebGenerator::isValidEdge(Graph& g, node u, node v) {
	auto isValid([&](node u, node v, float squaredDistance) {
		return ((squaredDistance <= neighRad * neighRad)
				&& (g.degree(u) <= maxNeigh)
				&& (g.degree(v) <= maxNeigh));
	});

	float x1 = g.getCoordinate(u, 0);
	float y1 = g.getCoordinate(u, 1);
	float x2 = g.getCoordinate(v, 0);
	float y2 = g.getCoordinate(v, 1);
	float sqrDist = squaredDistanceInUnitTorus(x1, y1, x2, y2);

	return isValid(u, v, sqrDist);
}

void DynamicPubWebGenerator::determineNeighborsOf(Graph& g, node u) {
	g.forNodes([&](node v) {
		if (isValidEdge(g, u, v)) {
			g.addEdge(u, v);
		}
	});
}

void DynamicPubWebGenerator::determineNeighbors(Graph& g) {
	g.forNodePairs([&](node u, node v) { // TODO: improve quadratic loop!
		if (isValidEdge(g, u, v)) {
			g.addEdge(u, v);
		}
	});
}

void DynamicPubWebGenerator::addNodesToArea(index area, count num, Graph& g) {
	Aux::RandomProbability randGen;

	for (index j = 0; j < num; ++j) {
		// compute random angle between [0, 2pi) and distance between [0, width/2]
		float angle = randGen.randomFloat() * 2.0 * PI;
		float dist = randGen.randomFloat() * denseAreaXYR[area].rad;

		// compute coordinates and adjust them
		float x = denseAreaXYR[area].x + cosf(angle) * dist;
		float y = denseAreaXYR[area].y + sinf(angle) * dist;
		moveNodeIntoUnitSquare(x, y);

		// create vertex with these coordinates
		g.addNode(x, y);
	}
}

void DynamicPubWebGenerator::fillDenseAreas(Graph& g) {
	Aux::RandomProbability randGen;

	for (index area = 0; area < numDenseAreas; ++area) {
		// choose center randomly, ensure complete cluster is within (0,1) without modifications
		denseAreaXYR[area].x = randGen.randomFloat();
		denseAreaXYR[area].y = randGen.randomFloat();
		addNodesToArea(area, numPerArea[area], g);
	}
}

void DynamicPubWebGenerator::chooseDenseAreaSizes() {
	denseAreaXYR.reserve(numDenseAreas);
	Aux::RandomProbability randGen;

	for (index area = 0; area < numDenseAreas; ++area) {
		// anti-quadratic probability distribution
		float f = randGen.randomFloat() * MIN_MAX_DENSE_AREA_FACTOR + 1.0f;
		denseAreaXYR[area].rad = (MAX_DENSE_AREA_RADIUS * f * f)
				/ (MIN_MAX_DENSE_AREA_FACTOR * MIN_MAX_DENSE_AREA_FACTOR);
	}
}

// randomly spread remaining vertices over whole area
void DynamicPubWebGenerator::spreadRemainingNodes(Graph& g) {
	Aux::RandomProbability randGen;

	while (g.numberOfNodes() < n) {
		float x = randGen.randomFloat();
		float y = randGen.randomFloat();
		g.addNode(x, y);
	}
}

// compute number of nodes per cluster, each cluster has approx. same density
void DynamicPubWebGenerator::chooseClusterSizes() {
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

Graph DynamicPubWebGenerator::generate() {
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
void DynamicPubWebGenerator::removeRandomNode(Graph& g) {
	Aux::RandomInteger randInt;
	node u = randInt.generate(0, (n - 1));
	g.removeNode(u);
}


#endif


} /* namespace NetworKit */
