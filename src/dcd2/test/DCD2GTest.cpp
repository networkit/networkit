/*
 * DCD2GTest.cpp
 *
 *  Created on: 07.01.2014
 *      Author: cls
 */

#include "DCD2GTest.h"

#include "../DynPLP.h"
#include "../DynPLM.h"
#include "../DynamicCommunityDetection.h"
#include "../../dynamics/DGSStreamParser.h"
#include "../../dynamics/GraphUpdater.h"

namespace NetworKit {


TEST_F(DCD2GTest, testGraphUpdater) {
	DGSStreamParser parser("input/example2.dgs", true, 1);
	DEBUG("getting event stream");
	std::vector<GraphEvent> stream = parser.getStream();
	Graph G;
	GraphUpdater gu(G);
	DEBUG("updating graph");
	gu.update(stream);

	EXPECT_EQ(3, G.numberOfNodes());

}


TEST_F(DCD2GTest, testDynPLP) {
	DGSStreamParser parser("input/example2.dgs");
	std::vector<GraphEvent> stream = parser.getStream();
	Graph G;
	DynPLP dynPLP;
	dynPLP.attachGraph(G);

	GraphUpdater gu(G);
	gu.update(stream);

	dynPLP.update(stream);
	Clustering zeta = dynPLP.detect();

	EXPECT_TRUE(zeta.isProper(G));

}

TEST_F(DCD2GTest, testDynPLPPrepStrategies) {
	DGSStreamParser parser("input/example2.dgs");
	std::vector<GraphEvent> stream = parser.getStream();

	Graph G;
	DynPLP dynPLP("isolate", 0);
	dynPLP.attachGraph(G);

	GraphUpdater gu(G);
	gu.update(stream);

	dynPLP.update(stream);
	Clustering zeta = dynPLP.detect();

	EXPECT_TRUE(zeta.isProper(G));

	Graph G2;
	DynPLP dynPLP2("isolateNeighbors", 0);
	dynPLP2.attachGraph(G2);

	GraphUpdater gu2(G2);
	gu2.update(stream);

	dynPLP2.update(stream);
	Clustering zeta2 = dynPLP2.detect();

	EXPECT_TRUE(zeta2.isProper(G));

}


TEST_F(DCD2GTest, testDynPLM) {
	DGSStreamParser parser("input/example2.dgs");
	std::vector<GraphEvent> stream = parser.getStream();
	Graph G;
	DynPLM dynPLM;
	dynPLM.attachGraph(G);

	GraphUpdater gu(G);
	gu.update(stream);

	dynPLM.update(stream);
	Clustering zeta = dynPLM.detect();

	EXPECT_TRUE(zeta.isProper(G));
}


TEST_F(DCD2GTest, testDynPLMPrepStrategies) {
	DGSStreamParser parser("input/example2.dgs");
	std::vector<GraphEvent> stream = parser.getStream();

	Graph G;
	DynPLM dynPLM("isolate", 0);
	dynPLM.attachGraph(G);

	GraphUpdater gu(G);
	gu.update(stream);

	dynPLM.update(stream);
	Clustering zeta = dynPLM.detect();

	EXPECT_TRUE(zeta.isProper(G));

	Graph G2;
	DynPLM dynPLM2("isolateNeighbors", 0);
	dynPLM2.attachGraph(G2);

	GraphUpdater gu2(G2);
	gu2.update(stream);

	dynPLM2.update(stream);
	Clustering zeta2 = dynPLM2.detect();

	EXPECT_TRUE(zeta2.isProper(G));

}


TEST_F(DCD2GTest, testDynamicCommunityDetectionWithPLP) {
	std::string path = "input/arxiv-qfin-author.dgs";
	DynamicCommunityDetection dynCD(path, "DynPLP", "isolate", 100);
	dynCD.run();

	INFO("quality timeline: " , Aux::vectorToString(dynCD.getTimeline("quality")));
}



TEST_F(DCD2GTest, testDynamicCommunityDetectionWithPLM) {
	std::string path = "input/arxiv-qfin-author.dgs";
	DynamicCommunityDetection dynCD(path, "DynPLM", "isolate", 100);
	dynCD.run();

	INFO("quality timeline: " , Aux::vectorToString(dynCD.getTimeline("quality")));

}


TEST_F(DCD2GTest, testDynamicCommunityDetectionWithStatic) {
	std::string path = "input/arxiv-qfin-author.dgs";
	DynamicCommunityDetection dynCD(path, "PLP", "isolate", 100);
	dynCD.run();

	INFO("quality timeline: " , Aux::vectorToString(dynCD.getTimeline("quality")));
}


TEST_F(DCD2GTest, tryDynPLMisolateNeighborsOnRealGraph) {
	std::string path = "input/arxiv-qfin-author.dgs";
	DynamicCommunityDetection dynCD(path, "DynPLM", "isolateNeighbors", 100, {"quality"});
	dynCD.run();

	INFO("quality timeline: " , Aux::vectorToString(dynCD.getTimeline("quality")));
}





} /* namespace NetworKit */
