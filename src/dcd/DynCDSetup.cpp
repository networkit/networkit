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
		nDetectors(dynDetectors.size()),
		tMax(tMax),
		deltaT(deltaT),
		staticAlgo(NULL),
		dynamicClusteringTimelines(nDetectors),
		qualityTimelines(nDetectors),
		nCommunitiesTimelines(nDetectors),
		continuityTimelines(nDetectors),
		checkMod(false),
		checkNumCom(false),
		checkSampledRand(false),
		checkNMID(false) {
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
	DynamicNMIDistance NMID;
	SampledRandMeasure sampledRand(500);

	// initialize graph
	gen->initializeGraph();


	auto runIt = [&](){
		// run algorithms and evaluation
		// inspect the current graph
		INFO("current graph : " << G->toString());

		nTimeline.push_back(G->numberOfNodes());
		mTimeline.push_back(G->numberOfEdges());

		// run the dynamic community detectors
		for (index detectorIndex = 0; detectorIndex < this->detectors.size(); ++detectorIndex) {
			DynamicCommunityDetector* dynCD = this->detectors.at(detectorIndex);

			INFO("running dynamic community detector " << dynCD->toString());
			dynamicClusteringTimelines.at(detectorIndex).push_back(dynCD->run());

			// TRACE("clustering looks like: " << Aux::vectorToString(results.at(detectorIndex).back().getVector()));
			assert (dynamicClusteringTimelines.at(detectorIndex).back().isProper(*G));

			// evaluations which need the current graph
			if (checkNumCom) {
				INFO("calculating number of clusters");
				count nCom = dynamicClusteringTimelines.at(detectorIndex).back().numberOfClusters();
				INFO("[RESULT] number of communities \t "<< detectors.at(detectorIndex)->toString() << "\t" << nCom);
				nCommunitiesTimelines.at(detectorIndex).push_back(nCom);
			}

			// modularity
			if (checkMod) {
				double mod = modularity.getQuality(dynamicClusteringTimelines.at(detectorIndex).back(), *G);
				INFO("[RESULT] modularity \t " << detectors.at(detectorIndex)->toString() << " \t " << mod);
				qualityTimelines.at(detectorIndex).push_back(mod);
			}

			// continuity by sampling
			if (checkSampledRand  && (dynamicClusteringTimelines[detectorIndex].size() >= 2)) {
				double cont = sampledRand.getDissimilarity(*G, dynamicClusteringTimelines.at(detectorIndex).at(dynamicClusteringTimelines.at(detectorIndex).size() - 2), dynamicClusteringTimelines.at(detectorIndex).back());
				INFO("[RESULT] continuity \t " << detectors.at(detectorIndex)->toString() << " \t " << cont);
				continuityTimelines.at(detectorIndex).push_back(cont);
			} else if (checkNMID && (dynamicClusteringTimelines.at(detectorIndex).size() >= 2)) {
				INFO("calculating continuity with NMID");
				Aux::Timer nmidTimer;
				nmidTimer.start();
				double cont = NMID.getDissimilarity(*G, dynamicClusteringTimelines.at(detectorIndex).at(dynamicClusteringTimelines.at(detectorIndex).size() - 2), dynamicClusteringTimelines.at(detectorIndex).back());
				nmidTimer.stop();
				INFO("calculating NMID took " << nmidTimer.elapsedTag());
				INFO("[RESULT] continuity NMID \t " << detectors.at(detectorIndex)->toString() << " \t " << cont);
				continuityTimelines.at(detectorIndex).push_back(cont);
			}
		} // end for multiple detectors

		// optionally also run a static community detector
		if (staticAlgo != NULL) {
			Aux::Timer staticRuntime;
			staticRuntime.start();
			//
			staticClusteringTimeline.push_back(staticAlgo->run(*G));
			//
			staticRuntime.stop();
			staticTimerTimeline.push_back(staticRuntime.elapsed().count());

			assert (staticClusteringTimeline.back().isProper(*G));

			if (checkNumCom) {

				INFO("calculating number of communities");
				count nCom = staticClusteringTimeline.back().numberOfClusters();
				INFO("[RESULT] number of communities \t " << staticAlgo->toString() << "\t " << nCom);
				this->staticNCommunitiesTimeline.push_back(nCom);
			}

			// modularity
			if (checkMod) {
				double mod = modularity.getQuality(staticClusteringTimeline.back(), *G);
				INFO("[RESULT] modularity \t" << staticAlgo->toString() << "\t " << mod);
				this->staticQualityTimeline.push_back(mod);
			}

			if (checkSampledRand  && (staticClusteringTimeline.size() >= 2)) { // if static algo has been set,
				double cont = sampledRand.getDissimilarity(*G, staticClusteringTimeline.at(staticClusteringTimeline.size() - 2), staticClusteringTimeline.back());
				INFO("[RESULT] continuity \t " << staticAlgo->toString() <<  " \t " << cont);
				this->staticContinuityTimeline.push_back(cont);
			} else if (checkNMID && (staticClusteringTimeline.size() >= 2)) {
				double cont = NMID.getDissimilarity(*G, staticClusteringTimeline.at(staticClusteringTimeline.size() - 2), staticClusteringTimeline.back());
				INFO("[RESULT] continuity NMID \t " << staticAlgo->toString() <<  " \t " << cont);
				this->staticContinuityTimeline.push_back(cont);
			}
		}

	};


	// for all community detectors, perform run
	bool sourceEnd = false;
	while (G->time() < tMax) {
		try {
			 // try generating the next batch of events
			INFO("receiving next batch of events");
			gen->generateTimeSteps(G->time() + deltaT);

			INFO("=========================== current time step: " << G->time() << " of " << tMax << " ========================================");

		} catch (std::logic_error& e) {
			INFO("exception caught: " << e.what());
			sourceEnd = true;
		}

		// RUN ALGORITHMS AND EVALUTION
		runIt();


		// break from the loop if source cannot produce more events
		if (sourceEnd) {
			INFO("breaking from loop because of end of source");
			break;
		}

	} // end while

	runtime.stop();

	INFO("=============================== RESULTS ===================================");

	INFO("setup runtime: " << runtime.elapsedTag());


	INFO("timeline \t n \t: " << Aux::vectorToString(this->nTimeline));
	INFO("timeline \t m \t: " << Aux::vectorToString(this->mTimeline));

	for (index detectorIndex = 0; detectorIndex < nDetectors; ++detectorIndex) {
		INFO("timeline \t " << detectors.at(detectorIndex)->toString() << " \t running time \t: " << Aux::vectorToString(this->detectors.at(detectorIndex)->getTimerHistory()));
		INFO("timeline \t " << detectors.at(detectorIndex)->toString() << " \t quality \t: " << Aux::vectorToString(this->qualityTimelines.at(detectorIndex)));
		INFO("timeline \t " << detectors.at(detectorIndex)->toString() << " \t # communities \t: " << Aux::vectorToString(this->nCommunitiesTimelines.at(detectorIndex)));
		INFO("timeline \t " << detectors.at(detectorIndex)->toString() << " \t continuity \t: " << Aux::vectorToString(this->continuityTimelines.at(detectorIndex)));
	}

	// static
	if (staticAlgo != NULL) {
		INFO("timeline \t " << this->staticAlgo->toString() << " \t running time \t: " << Aux::vectorToString(this->staticTimerTimeline));
		INFO("timeline \t " << this->staticAlgo->toString() << " \t quality \t: " << Aux::vectorToString(this->staticQualityTimeline));
		INFO("timeline \t " << this->staticAlgo->toString() << " \t # communities \t: " << Aux::vectorToString(this->staticNCommunitiesTimeline));
		INFO("timeline \t " << this->staticAlgo->toString() << " \t continuity \t: " << Aux::vectorToString(this->staticContinuityTimeline));
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

void DynCDSetup::checkContinuity() {
	// this->checkSampledRand = true;
	this->checkNMID = true;
}

} /* namespace NetworKit */
