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

void DynamicCommunityDetection::run(std::string inputPath, std::string algoName) {
	DGSStreamParser parser(inputPath, true, 1);
	DEBUG("getting event stream");
	std::vector<GraphEvent> stream = parser.getStream();
	Graph G;
	GraphUpdater gu(G);
	DEBUG("updating graph");
	gu.update(stream);

	DynCommunityDetector* algo = NULL;
	if (algoName == "DynPLP") {
		algo = new DynPLP();
	} else if (algoName == "DynPLM") {
		algo = new DynPLM();
	}

	algo->attachGraph(G);


}

} /* namespace NetworKit */
