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
		ERROR("deltaT >= tMax");
		throw std::runtime_error("deltaT must be smaller than tMax");
	} else {
		INFO("deltaT " << deltaT);
		INFO("tMax " << tMax);

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
	Aux::Timer runtime;
	runtime.start();

	// helpers
	Modularity modularity;

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
			INFO("G->time() < tMax is TRUE");

			gen->generateTimeSteps(G->time() + deltaT);
			// inspect the current graph
			INFO("current graph : " << G->toString());

			// run the dynamic community detectors
			for (count i = 0; i < this->detectors.size(); ++i) {
				DynamicCommunityDetector* dynCD = this->detectors[i];
				INFO("running dynamic community detector " << dynCD->toString());
				results[i].push_back(dynCD->run());

				// evaluations which need the current graph
				if (checkMod) {
					double mod = modularity.getQuality(results[i].back(), *G);
					INFO("found communities have modularity: " << mod);
				}
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

	} // end while

	runtime.stop();
	INFO("setup runtime: " << runtime.elapsedTag());

	if (checkNumCom) {
		for (std::vector<Clustering> dynZeta : results) {
			for (Clustering zeta : dynZeta) {
				INFO("calculating number of clusters");
				DEBUG("number of clusters: " << zeta.numberOfClusters());
			}
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

void DynCDSetup::checkModularity() {
	this->checkMod = true;
}

void DynCDSetup::checkNumberOfCommunities() {
	this->checkNumCom = true;
}

} /* namespace NetworKit */
