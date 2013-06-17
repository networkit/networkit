/*
 * DynCDSetup.cpp
 *
 *  Created on: 16.06.2013
 *      Author: cls
 */

#include "DynCDSetup.h"

namespace NetworKit {

DynCDSetup::DynCDSetup(DynamicGraphGenerator& dynGen, std::vector<DynamicCommunityDetector*>& dynDetectors, count tMax, count deltaT) :
		gen(&dynGen),
		detectors(dynDetectors),
		tMax(tMax),
		deltaT(deltaT) {
	if (deltaT >= tMax) {
		throw std::runtime_error("deltaT must be smaller than tMax");
	}

	// create graph and proxy instances
	Gproxy = this->gen->newGraph();
	G = Gproxy->G;

	DEBUG("setting up " << detectors.size() << " community detection algorithms");
	// register algorithms
	for (DynamicCommunityDetector* dynCD : detectors) {
		Gproxy->registerObserver(dynCD); 	// register community detection algorithm as observer
		dynCD->setGraph(*G);				// point the community detection algorithm to the graph
	}

}

DynCDSetup::~DynCDSetup() {
	// TODO Auto-generated destructor stub
}

void DynCDSetup::run() {

	// initialize graph
	gen->initializeGraph();

	// store the resulting clusterings in here
	std::vector<std::vector<Clustering> > results;
	for (count i = 0; i < this->detectors.size(); ++i) {
		std::vector<Clustering> dynZeta;
		results.push_back(dynZeta);
	}


	// for all community detectors, perform run

	while (G->time() <= tMax) {
		gen->generateTimeSteps(G->time() + deltaT);
		if (G->time() % deltaT == 0) {
			for (count i = 0; i < this->detectors.size(); ++i) {
				DynamicCommunityDetector* dynCD = this->detectors[i];
				INFO("running dynamic community detector " << dynCD->toString());
				results[i].push_back(dynCD->run());
			}
		}

	}
	for (std::vector<Clustering> dynZeta : results) {
		for (Clustering zeta : dynZeta) {
			DEBUG("number of clusters: " << zeta.numberOfClusters());
		}
	}

}

} /* namespace NetworKit */
