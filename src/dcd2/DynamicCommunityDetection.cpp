/*
 * DynamicCommunityDetection.cpp
 *
 *  Created on: 10.01.2014
 *      Author: cls
 */

#include "DynamicCommunityDetection.h"

#include "DynPLP.h"
#include "DynPLM.h"
#include "StaticAdapter.h"
#include "../dynamics/DGSStreamParser.h"
#include "../dynamics/GraphUpdater.h"
#include "../clustering/Modularity.h"
#include "../auxiliary/Timer.h"
#include "../community/PLP.h"
#include "../community/PLM2.h"

namespace NetworKit {

DynamicCommunityDetection::DynamicCommunityDetection(std::string inputPath, std::string algoName, count interval, bool recordQuality) : inputPath(inputPath), algoName(algoName), interval(interval), recordQuality(recordQuality) {
}

void DynamicCommunityDetection::run() {
	DEBUG("reading full event stream from file");
	DGSStreamParser parser(inputPath, true, 1);
	std::vector<GraphEvent> stream = parser.getStream();
	DEBUG("setting up graph");
	Graph G;
	GraphUpdater gu(G);

	DEBUG("setting up algorithms");
	DynCommunityDetector* algo = NULL;
	if (algoName == "DynPLP") {
		algo = new DynPLP();
	} else if (algoName == "DynPLM") {
		algo = new DynPLM();
	} else if (algoName == "PLP") {
		algo = new StaticAdapter(new PLP());
	} else if (algoName == "PLM") {
		algo = new StaticAdapter(new PLM2());
	}

	DEBUG("attaching graph");
	algo->attachGraph(G);

	Aux::Timer timer;


	// slicing stream by timesteps
	auto i = stream.begin();
	auto j = stream.begin();

	// advance the iterators forward to the next <interval> time steps
	auto advance = [&](){
		count steps = 0;
		i = j;
		for (; (steps < interval) && (j < stream.end()); ++j) {
			if (j->type == GraphEvent::TIME_STEP) {
				steps++;
			}
		}
	};

	DEBUG("advancing dynamic graph");
	while (j != stream.end()) {
		advance(); // advance iterators
		std::vector<GraphEvent> slice(i, j);
		DEBUG("updating graph");
		gu.update(slice);
		DEBUG("number of nodes: " << G.numberOfNodes());


		timer.start();
		//
		algo->update(slice);
		//
		timer.stop();
		updateTime.push_back(timer.elapsedMilliseconds());

		timer.start();
		//
		Clustering zeta = algo->detect();
		//
		timer.stop();
		detectTime.push_back(timer.elapsedMilliseconds());

		DEBUG("upper community id bound: " << zeta.upperBound());
		// DEBUG("zeta at time " << G.time() << ": " << Aux::vectorToString(zeta.getVector()));

		if (recordQuality) {
			Modularity mod;
			quality.push_back(mod.getQuality(zeta, G));
			INFO("modularity at time " << G.time() << ": " << quality.back());
		}
	}


	// record graph size
	size = gu.getSizeTimeline();

}

std::vector<count> DynamicCommunityDetection::getUpdateTimeline() {
	return updateTime;
}

std::vector<count> DynamicCommunityDetection::getDetectTimeline() {
	return detectTime;
}

std::vector<double> DynamicCommunityDetection::getQualityTimeline() {
	return quality;
}

std::vector<double> DynamicCommunityDetection::getContinuityTimeline() {
	return continuity;
}

std::vector<std::pair<count, count> > DynamicCommunityDetection::getGraphSizeTimeline() {
	return size;
}

} /* namespace NetworKit */
