/*
 * InputGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "InputGTest.h"

namespace EnsembleClustering {

TEST_F(InputGTest, testGraphIOEdgeList) {
	GraphGenerator graphGen;
	Graph G = graphGen.makeCircularGraph(20);
	GraphIO graphio;
	std::string path = "sandbox/edgelist.txt";
	graphio.writeEdgeList(G, path);

	bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(InputGTest, testGraphIOAdjacencyList) {
	GraphGenerator graphGen;
	Graph G = graphGen.makeCircularGraph(20);
	GraphIO graphio;
	std::string path = "sandbox/circular.adjlist";
	graphio.writeAdjacencyList(G, path);

	bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "A file should have been created : " << path;
}


TEST_F(InputGTest, testGraphIOForIsolatedNodes) {
	GraphGenerator graphGen;
	Graph G(20);
	GraphIO graphio;
	std::string path = "sandbox/isolated.adjlist";
	graphio.writeAdjacencyList(G, path);

	bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "A file should have been created : " << path;
}



TEST_F(InputGTest, testMETISGraphReader) {
	std::string path = "input/jazz.graph";

	METISGraphReader reader;
	Graph G = reader.read(path);

	EXPECT_FALSE(G.isEmpty());
}




TEST_F(InputGTest, testClusteringWriterAndReader) {
	// write clustering first
	std::string path = "sandbox/example.clust";

	GraphGenerator graphGen;
	count n = 100;
	count k = 3;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Clustering zeta = clusteringGen.makeRandomClustering(G, k);

	ClusteringWriter writer;
	writer.write(zeta, path);

	// check if file exists
	bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "clustering file should have been written to: " << path;


	ClusteringReader reader;
	Clustering read = reader.read(path);

	EXPECT_EQ(n, read.numberOfNodes()) << "read clustering should contain n nodes";
	EXPECT_TRUE(read.isProper(G)) << "read clustering should be proper clustering of G";
	EXPECT_TRUE(read.equals(zeta, G)) << "read clustering should be identical to created clustering";
}



} /* namespace EnsembleClustering */
