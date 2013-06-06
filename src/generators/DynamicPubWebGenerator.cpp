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
		DynamicGraphGenerator(proxy),
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
	DEBUG("G: " << this->Gproxy << ", number of nodes: " << Gproxy->G->numberOfNodes());
//	delete Gproxy->G;
	Gproxy->G = new Graph(staticGen.generate());
	G = Gproxy->G;
	DEBUG("Gproxy: " << this->Gproxy << ", number of nodes: " << Gproxy->G->numberOfNodes());
	DEBUG("G:      " << this->G << ", number of nodes: " << G->numberOfNodes());
}

node DynamicPubWebGenerator::addNode() {

	// identify dense area (or remaining)
	// -> find out where rand value falls into the interval divided by
	//    the fraction of the cluster size, don't forget the non-dense area
	// TODO
	std::vector<float> probInterval(staticGen.numDenseAreas+1);
	count n = Gproxy->G->numberOfNodes();

	for (index i = 0; i < staticGen.numDenseAreas; ++i) {
		probInterval[i] = (float) staticGen.numPerArea[i] / (float) n;
	}
	// prefix sums
	for (index i = 1; i < staticGen.numDenseAreas; ++i) {
		probInterval[i] += probInterval[i-1];
	}
	probInterval[staticGen.numDenseAreas] = 1.0;

	Aux::RandomProbability randProb;
	float r = randProb.randomFloat();

	index area = 0;
	while (r > probInterval[area] && area < staticGen.numDenseAreas) {
		++area;
	}

	// identify random location in that area
	count dims = 2;
	count num = 1;
	staticGen.addNodesToArea(area, num, *G);

	// insert edges according to rules in determineNeighbors
	node u = n; // inserted last
	staticGen.determineNeighborsOf(*G, u);

	return u;
}


void DynamicPubWebGenerator::generateWhile(std::function<bool(void)> cont) {
	// TODO: delete, so far only add

	// incoming action: add or delete vertex
	// => add or delete edges accordingly

	while (cont()) {
		node u = this->addNode();
		TRACE("adding node " << u);
	}

}

} /* namespace NetworKit */
