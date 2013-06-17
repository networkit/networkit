/*
 * DynamicPubWebGenerator.cpp
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#include "DynamicPubWebGenerator.h"

namespace NetworKit {

DynamicPubWebGenerator::DynamicPubWebGenerator(count numInitialNodes, count numberOfDenseAreas, float neighborhoodRadius,
		count maxNumberOfNeighbors) :
		DynamicGraphSource(),
		staticGen(numInitialNodes, numberOfDenseAreas, neighborhoodRadius,
				maxNumberOfNeighbors)
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
	staticGen.chooseDenseAreaSizes();
	staticGen.chooseClusterSizes();
	staticGen.fillDenseAreas(*G);
	staticGen.spreadRemainingNodes(*G);
	staticGen.determineNeighbors(*G);
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
		this->staticGen.moveNodeIntoUnitSquare(x, y);

		// move node
		this->G->setCoordinate(u, 0, x);
		this->G->setCoordinate(u, 1, y);
	});
}


void DynamicPubWebGenerator::generate() {
	// nodes move to new position
	moveNodesRandomly();

	// decide if new edges have to be inserted
	// FIXME: inefficient!
	G->forNodePairs([&](node u, node v) {
		bool valid = staticGen.isValidEdge(*G, u, v);
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
}

} /* namespace NetworKit */
