/*
 * PropertiesGTest_Ritter.cpp
 */

#ifndef NOGTEST

#include "PropertiesGTest_Ritter.h"

#include "../CoreDecomposition_Ritter.h"
#include "../CoreDecomposition_Ritter_ownList.h"


namespace NetworKit {

PropertiesGTest_Ritter::PropertiesGTest_Ritter() {

}

PropertiesGTest_Ritter::~PropertiesGTest_Ritter() {

}

TEST_F(PropertiesGTest_Ritter, testCoreDecomposition) {
	count n = 16;
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

TEST_F(PropertiesGTest_Ritter, testCoreDecompositionWithOwnListImpl) {
	count n = 16;
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
	CoreDecomposition_Ritter_ownList coreDec;
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
