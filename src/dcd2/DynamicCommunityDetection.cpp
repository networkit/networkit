/*
 * DynamicCommunityDetection.cpp
 *
 *  Created on: 10.01.2014
 *      Author: cls
 */

#include "DynamicCommunityDetection.h"

#include "DynPLP.h"
#include "DynPLM.h"
#include "../dynamics/DGSStreamParser.h"
#include "../dynamics/GraphUpdater.h"
#include "../clustering/Modularity.h"

namespace NetworKit {

DynamicCommunityDetection::DynamicCommunityDetection() {
	// TODO Auto-generated constructor stub

}

void DynamicCommunityDetection::run(std::string inputPath, std::string algoName, count interval) {
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
	}

	algo->attachGraph(G);

	// slicing stream by timesteps
	auto i = stream.begin();
	auto j = stream.begin();


	auto advance = [&](){
		count steps = 0;
		i = j;
		for (; (steps < interval) && (j < stream.end()); ++j) {
			if (j->type == GraphEvent::TIME_STEP) {
				steps++;
			}
		}
	};



	while (j != stream.end()) {
		advance(); // advance iterators
		std::vector<GraphEvent> slice(i, j);
		DEBUG("updating graph");
		gu.update(slice);
		DEBUG("number of nodes: " << G.numberOfNodes());

		algo->process(slice);
		Clustering zeta = algo->retrieve();
		DEBUG("upper community id bound: " << zeta.upperBound());
		// DEBUG("zeta at time " << G.time() << ": " << Aux::vectorToString(zeta.getVector()));

		Modularity mod;
		INFO("modularity at time " << G.time() << ": " << mod.getQuality(zeta, G));

	}




}

} /* namespace NetworKit */
