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

	// slicing stream
	auto i = stream.begin();
	auto j = i + interval;
	if (j > stream.end()) {
		j = stream.end();
	}

	while (i < stream.end()) { // TODO: rest of stream
		std::vector<GraphEvent> slice(i, j);
		DEBUG("updating graph");
		gu.update(slice);
		DEBUG("number of nodes: " << G.numberOfNodes());

		algo->process(slice);
		Clustering zeta = algo->retrieve();


		//
		i = j;
		j = i + interval;
		if (j > stream.end()) {
			j = stream.end();
		}
	}




}

} /* namespace NetworKit */
