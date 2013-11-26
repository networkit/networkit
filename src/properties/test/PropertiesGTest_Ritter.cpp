/*
 * PropertiesGTest_Ritter.cpp
 */

#ifndef NOGTEST

#include "PropertiesGTest_Ritter.h"

#include "../CoreDecomposition_Ritter.h"
#include "../GraphProperties_Ritter.h"

namespace NetworKit {

PropertiesGTest_Ritter::PropertiesGTest_Ritter() {

}

PropertiesGTest_Ritter::~PropertiesGTest_Ritter() {

}

bool goodDiameterEstimate(count diameter, count lowerBound, count upperBound) {
	return lowerBound <= diameter &&
			upperBound >= diameter &&
			(diameter < 3 || lowerBound > diameter - 3) &&
			upperBound < diameter + 3;
}

TEST_F(PropertiesGTest_Ritter, testEstimateDiameterRangeEmptyGraph) {
	Graph G(2);
	EXPECT_EQ(0, GraphProperties::estimatedDiameterRange_Ritter(G).first);
	EXPECT_EQ(0, GraphProperties::estimatedDiameterRange_Ritter(G).second);
}

TEST_F(PropertiesGTest_Ritter, testEstimateDiameterRangeGraphNotConnected) {
	Graph G(3);
	G.addEdge(0, 1);
	EXPECT_EQ(GraphProperties_Ritter::INF_DIST, GraphProperties::estimatedDiameterRange_Ritter(G).first);
	EXPECT_EQ(GraphProperties_Ritter::INF_DIST, GraphProperties::estimatedDiameterRange_Ritter(G).second);

}

TEST_F(PropertiesGTest_Ritter, testEstimateDiameterRange) {
	Graph G(2);
	std::pair<count, count> bounds;

	G.addEdge(0, 1);
	bounds = GraphProperties::estimatedDiameterRange_Ritter(G);
	EXPECT_PRED3(goodDiameterEstimate, 1, bounds.first, bounds.second) << "Graph with 1 edge";

	G.addNode();
	G.addEdge(1, 2);
	bounds = GraphProperties::estimatedDiameterRange_Ritter(G);
	EXPECT_PRED3(goodDiameterEstimate, 2, bounds.first, bounds.second) << "Graph with 2 edges";

	G.addEdge(0, 2);
	bounds = GraphProperties::estimatedDiameterRange_Ritter(G);
	EXPECT_PRED3(goodDiameterEstimate, 1, bounds.first, bounds.second) << "Graph with 3 edges";

	G.addNode();
	G.addEdge(0, 3);
	bounds = GraphProperties::estimatedDiameterRange_Ritter(G);
	EXPECT_PRED3(goodDiameterEstimate, 2, bounds.first, bounds.second) << "Graph with 4 edges";

	G.addNode();
	G.addEdge(3, 4);
	bounds = GraphProperties::estimatedDiameterRange_Ritter(G);
	EXPECT_PRED3(goodDiameterEstimate, 3, bounds.first, bounds.second) << "Graph with 5 edges";

	G.addNode();
	G.addEdge(4, 5);
	bounds = GraphProperties::estimatedDiameterRange_Ritter(G);
	EXPECT_PRED3(goodDiameterEstimate, 4, bounds.first, bounds.second) << "Graph with 6 edges";
}

TEST_F(PropertiesGTest_Ritter, testDiameter) {
	Graph G(0);

	EXPECT_EQ(0, GraphProperties_Ritter::diameter(G)) << "Graph with 0 edges";

	G.addNode();
	G.addNode();
	G.addEdge(0, 1);
	EXPECT_EQ(1, GraphProperties_Ritter::diameter(G)) << "Graph with 1 edge";

	G.addEdge(1, 2);
	EXPECT_EQ(2, GraphProperties_Ritter::diameter(G)) << "Graph with 2 edges";

	G.addEdge(0, 2);
	EXPECT_EQ(1, GraphProperties_Ritter::diameter(G)) << "Graph with 3 edges";

	G.addNode();
	G.addEdge(0, 3);
	EXPECT_EQ(2, GraphProperties_Ritter::diameter(G)) << "Graph with 4 edges";

	G.addNode();
	G.addEdge(3, 4);
	EXPECT_EQ(3, GraphProperties_Ritter::diameter(G)) << "Graph with 5 edges";

	G.addNode();
	G.addEdge(4, 5);
	EXPECT_EQ(4, GraphProperties_Ritter::diameter(G)) << "Graph with 6 edges";

	G.addNode();
	G.addEdge(6, 7);
	EXPECT_EQ(4, GraphProperties_Ritter::diameter(G)) << "Graph with 7 edges";

	G.addNode();
	G.addNode();
	G.addNode();
	G.addEdge(7, 8);
	G.addEdge(8, 9);
	G.addEdge(9, 10);
	EXPECT_EQ(4, GraphProperties_Ritter::diameter(G)) << "Graph with 10 edges";
	
	G.addNode();
	G.addEdge(10, 11);
	EXPECT_EQ(5, GraphProperties_Ritter::diameter(G)) << "Graph with 11 edges";

	G.addEdge(8, 0);
	EXPECT_EQ(7, GraphProperties_Ritter::diameter(G)) << "Graph with 12 edges";
}

TEST_F(PropertiesGTest_Ritter, testCoreDecomposition) {
	const count n = 16;
	Graph G(n);

	// create graph used in Baur et al. and network analysis lecture
	G.addEdge(2, 4);
	G.addEdge(3, 4);
	G.addEdge(4, 5);
	G.addEdge(5, 7);
	G.addEdge(6, 7);

	G.addEdge(6, 8);
	G.addEdge(6, 9);
	G.addEdge(6, 11);
	G.addEdge(7, 12);
	G.addEdge(8, 9);

	G.addEdge(8, 10);
	G.addEdge(8, 11);
	G.addEdge(8, 13);
	G.addEdge(9, 10);
	G.addEdge(9, 11);

	G.addEdge(9, 13);
	G.addEdge(10, 11);
	G.addEdge(10, 13);
	G.addEdge(10, 14);
	G.addEdge(11, 13);

	G.addEdge(11, 14);
	G.addEdge(12, 15);
	G.addEdge(13, 14);
	G.addEdge(14, 15);

	EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
	EXPECT_EQ(24, G.numberOfEdges()) << "should have 24 edges";

	// compute core decomposition
	CoreDecomposition_Ritter coreDec;
	std::vector<count> coreness = coreDec.run(G);

	EXPECT_EQ(16, coreness.size()) << "coreness must be set for every node";
	EXPECT_EQ(0, coreness[0]) << "expected coreness";
	EXPECT_EQ(0, coreness[1]) << "expected coreness";
	EXPECT_EQ(1, coreness[2]) << "expected coreness";
	EXPECT_EQ(1, coreness[3]) << "expected coreness";
	EXPECT_EQ(1, coreness[4]) << "expected coreness";
	EXPECT_EQ(1, coreness[5]) << "expected coreness";
	EXPECT_EQ(3, coreness[6]) << "expected coreness";
	EXPECT_EQ(2, coreness[7]) << "expected coreness";
	EXPECT_EQ(4, coreness[8]) << "expected coreness";
	EXPECT_EQ(4, coreness[9]) << "expected coreness";
	EXPECT_EQ(4, coreness[10]) << "expected coreness";
	EXPECT_EQ(4, coreness[11]) << "expected coreness";
	EXPECT_EQ(2, coreness[12]) << "expected coreness";
	EXPECT_EQ(4, coreness[13]) << "expected coreness";
	EXPECT_EQ(3, coreness[14]) << "expected coreness";
	EXPECT_EQ(2, coreness[15]) << "expected coreness";
}

} /* namespace NetworKit */

#endif /*NOGTEST*/
