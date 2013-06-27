/*
 * DynCDSetup.cpp
 *
 *  Created on: 16.06.2013
 *      Author: cls
 */

#include "DynCDSetup.h"

namespace NetworKit {

DynCDSetup::DynCDSetup(DynamicGraphSource& dynGen, std::vector<DynamicCommunityDetector*>& dynDetectors, count tMax, count deltaT) :
		gen(&dynGen),
		detectors(dynDetectors),
		tMax(tMax),
		deltaT(deltaT),
		staticAlgo(NULL) {
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

	// store the resulting clusterings in results

	for (count i = 0; i < this->detectors.size(); ++i) {
		std::vector<Clustering> dynZeta;
		results.push_back(dynZeta);
	}


	// for all community detectors, perform run
	while (G->time() < tMax) {
		INFO("time: " << G->time() << " of " << tMax);
		try {
			gen->generateTimeSteps(G->time() + deltaT);
			// run the dynamic community detectors
			for (count i = 0; i < this->detectors.size(); ++i) {
				DynamicCommunityDetector* dynCD = this->detectors[i];
				INFO("running dynamic community detector " << dynCD->toString());
				results[i].push_back(dynCD->run());

			}
			// optionally also run a static community detector
			if (staticAlgo != NULL) {
				staticClusterings.push_back(staticAlgo->run(*G));
			}
		} catch (std::logic_error& e) {
			INFO("source cannot produce any more events");
			// perform algorithm runs for the last time
			for (count i = 0; i < this->detectors.size(); ++i) {
				DynamicCommunityDetector* dynCD = this->detectors[i];
				INFO("running dynamic community detector " << dynCD->toString());
				results[i].push_back(dynCD->run());
			}
			break;
		}

	}
	for (std::vector<Clustering> dynZeta : results) {
		for (Clustering zeta : dynZeta) {
			DEBUG("number of clusters: " << zeta.numberOfClusters());
		}
	}

	for (DynamicCommunityDetector* dynCD : this->detectors){
		INFO("timer history for algorithm " << dynCD->toString() << ": " << Aux::vectorToString(dynCD->getTimerHistory()));
	}

}

Graph* DynCDSetup::getGraph() {
	return this->G;
}

void DynCDSetup::setStatic(Clusterer* staticAlgo) {
	this->staticAlgo = staticAlgo;
}

Graph DynCDSetup::getGraphCopy() {
	return *(this->G);
}

} /* namespace NetworKit */
