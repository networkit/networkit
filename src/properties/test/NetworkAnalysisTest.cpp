/*
 * NetworkAnalysisTest.cpp
 *
 *  Created on: 18.11.2013
 *      Author: christianocker
 */

#ifndef NOGTEST

#include "NetworkAnalysisTest.h"


namespace NetworKit {

TEST_F(NetworkAnalysisTest, testCoreDecomposition) {
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
	CoreDecomposition coreDec;
	std::vector<count> coreness = coreDec.run(G);

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

TEST_F(NetworkAnalysisTest, tryCoreDecompositionDimacsGraphs) {
  METISGraphReader reader = METISGraphReader();

  std::vector<std::string> graphPaths;
  graphPaths.push_back("dimacs/celegans_metabolic.graph");
  graphPaths.push_back("dimacs/polblogs.graph");
  graphPaths.push_back("dimacs/hep-th.graph");

  for(std::string path: graphPaths) {
    std::cout << "Graph: " << path << std::endl;
    Graph G = reader.read(std::string("input/") + path);

    CoreDecomposition coreDec;
    std::vector<count> coreness = coreDec.run(G);
    
    std::ofstream solutionFile(std::string("output/") + path.substr(0, path.size() - 6) + ".sol", std::ofstream::out);
    for(count c: coreness) {
      solutionFile << c << std::endl;
    }
  }
}

TEST_F(NetworkAnalysisTest, testClusteringCoefficient) {
	count n = 10;
	Graph G(n);

  G.addEdge(0, 1);
  G.addEdge(1, 2);
  G.addEdge(2, 3);
  G.addEdge(3, 4);
  G.addEdge(4, 5);

	G.addEdge(6, 7);
	G.addEdge(6, 8);
	G.addEdge(6, 9);
	G.addEdge(7, 8);
	G.addEdge(7, 9);
  G.addEdge(8, 9);

	EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
	EXPECT_EQ(11, G.numberOfEdges()) << "should have 11 edges";

  GlobalClusteringCoefficient cc;
  EXPECT_NEAR(0.75, cc.approximate(G, 100000), 0.01);
}

TEST_F(NetworkAnalysisTest, tryClusteringCoefficientDimacsGraphs) {
  METISGraphReader reader = METISGraphReader();

  std::vector<std::string> graphPaths;
  graphPaths.push_back("input/dimacs/celegans_metabolic.graph");
  graphPaths.push_back("input/dimacs/polblogs.graph");
  graphPaths.push_back("input/dimacs/hep-th.graph");

  for(std::string path: graphPaths) {
    std::cout << "Graph: " << path << std::endl;
    Graph G = reader.read(path);

    GlobalClusteringCoefficient cc;
    std::cout << "Global cluster coefficient: " << cc.approximate(G, 10000000) << std::endl;
  }
}

} /* namespace NetworKit */

#endif /*NOGTEST*/
